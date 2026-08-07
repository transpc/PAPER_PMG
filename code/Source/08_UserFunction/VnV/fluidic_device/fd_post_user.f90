!
      SUBROUTINE fd_post_user
!
!     User-defined post-process for Fludic Device problem
!    
      USE VOL_DATA        , ONLY: cell
      USE Zmpi            , ONLY: maxmt_fluid
      USE Zzone           , ONLY: ncell_fluid
      USE Zcore           , ONLY: myrank,np 
      USE Zparam          , ONLY: ndim
      USE Znum_cell       , ONLY: i_neigh,neigh
      USE Zbc_index       , ONLY: nbcon
      USE Zconst2         , ONLY: dt
      USE Ztimecon        , ONLY: time
      USE Zcoord1         , ONLY: xloc
      USE Zcoord3         , ONLY: vol
      USE Zfluidic_device , ONLY: flow_si,flow_sp,flow_fd
      USE Zpress          , ONLY: p
      USE Zvec_major      , ONLY: flux_l_nf,flux_g_nf
      USE Zvec_geo        , ONLY: xn_nf
!
      IMPLICIT NONE
!
!.....External function
      INTEGER :: get_nf_number
!.....Local variables
      INTEGER :: i,j,k,j0
      INTEGER :: ix
      INTEGER :: nf_number,i1
      INTEGER, SAVE :: num_f1,num_f2,num_f3,n_top,n_top_sp,n_top_fd,n_out1,n_out2
      LOGICAL, SAVE :: INITIAL=.TRUE.
      REAL(8) :: volume,height,area,volume_sp,level_sp,area_sp
      REAL(8) :: p_top,t_top
      REAL(8) :: flux_l,flux_g
      REAL(8), SAVE :: time1      
!.....Local arrays
      REAL(8) :: tmp(3)
      REAL(8) :: xn(maxmt_fluid,3)
      INTEGER, SAVE :: nf1(2,20),nf2(3,30),nf3(2,20)
!
      IF(INITIAL)THEN
!   bug time1 not initialized
         time1=0.d0
         IF(myrank.eq.0)THEN
!           OPEN(999,file='fd.dat')
            OPEN(333,file='VD6_ref.dat')
         ENDIF
         CALL get_scalar_variable_n_ndim(xn_nf,xn)
         num_f1=0
         num_f2=0
         num_f3=0
         DO i=1,ncell_fluid
!
            IF(xloc(i,1).gt.1.04115.and.xloc(i,1).lt.1.38995.and.  &
               xloc(i,2).gt.1.04115.and.xloc(i,2).lt.1.38995)THEN
               IF(xloc(i,ndim).gt.0.125.and.xloc(i,ndim).lt.0.25)THEN
                  num_f1=num_f1+1
                  nf1(1,num_f1)=i
                   DO j=i_neigh(i),i_neigh(i+1)-1
                     IF(xn(j,3).gt.0.5d0) nf1(2,num_f1)=j-i_neigh(i)+1
                  ENDDO
               ENDIF
               IF(xloc(i,ndim).gt.0..and.xloc(i,ndim).lt.0.125)THEN
                  num_f3=num_f3+1
                  nf3(1,num_f3)=i
                  DO j=i_neigh(i),i_neigh(i+1)-1
                     IF(xn(j,3).lt.-0.5d0) nf3(2,num_f3)=j-i_neigh(i)+1
                  ENDDO
               ENDIF
!
               IF(xloc(i,ndim).gt.11.2d0) n_top=i
               IF(xloc(i,ndim).gt.4.58d0 .and.xloc(i,ndim).lt.4.81d0) n_top_sp=i
               IF(xloc(i,ndim).gt.0.125  .and.xloc(i,ndim).lt.0.25d0) n_out1=i
               IF(xloc(i,ndim).gt.-1.73d0.and.xloc(i,ndim).lt.-1.15d00) n_out2=i
            ENDIF
!
            IF(xloc(i,ndim).gt.0..and.xloc(i,ndim).lt.0.25)THEN
               IF(xloc(i,1).gt.0.899.and.xloc(i,1).lt.1.04.and.xloc(i,2).gt.1.04.and.xloc(i,2).lt.1.39)THEN
                  num_f2=num_f2+1
                  nf2(1,num_f2)=i
                   DO j=i_neigh(i),i_neigh(i+1)-1
                     IF(xn(j,1).gt.0.5d0) nf2(2,num_f2)=j-i_neigh(i)+1
                  ENDDO
               ENDIF
               IF(xloc(i,1).gt.1.39.and.xloc(i,1).lt.1.53.and.xloc(i,2).gt.1.04.and.xloc(i,2).lt.1.39)THEN
                  num_f2=num_f2+1
                  nf2(1,num_f2)=i
                  DO j=i_neigh(i),i_neigh(i+1)-1
                     IF(xn(j,1).lt.-0.5d0) nf2(2,num_f2)=j-i_neigh(i)+1
                  ENDDO
               ENDIF
               IF(xloc(i,2).gt.0.899.and.xloc(i,2).lt.1.04.and.xloc(i,1).gt.1.04.and.xloc(i,1).lt.1.39)THEN
                  num_f2=num_f2+1
                  nf2(1,num_f2)=i
                  DO j=i_neigh(i),i_neigh(i+1)-1
                     IF(xn(j,2).gt.0.5d0) nf2(2,num_f2)=j-i_neigh(i)+1
                  ENDDO
               ENDIF
               IF(xloc(i,2).gt.1.39.and.xloc(i,2).lt.1.53.and.xloc(i,1).gt.1.04.and.xloc(i,1).lt.1.39)THEN
                  num_f2=num_f2+1
                  nf2(1,num_f2)=i
                  DO j=i_neigh(i),i_neigh(i+1)-1
                     IF(xn(j,2).lt.-0.5d0) nf2(2,num_f2)=j-i_neigh(i)+1
                  ENDDO
               ENDIF
            ENDIF
!
            IF(xloc(i,1).gt.0.444.and.xloc(i,1).lt.0.616.and.  &
               xloc(i,2).gt.1.13.and.xloc(i,2).lt.1.3)THEN
               IF(xloc(i,ndim).gt.0.636.and.xloc(i,ndim).lt.0.793) n_top_fd=i
            ENDIF
         ENDDO
         INITIAL=.FALSE.
      ENDIF
!
      flow_sp=0.d0
      flow_fd=0.d0
      flow_si=0.d0
      DO ix=1,num_f1
         i=nf1(1,ix)
         j0=i_neigh(i)-1
         j=nf1(2,ix)
         nf_number=get_nf_number(nbcon(j+j0))
         IF(nf_number.le.3) THEN
            CALL get_vector_disp(j,i,i1)
            IF(i1.gt.0) THEN
               flux_l=flux_l_nf(i1)
               flux_g=flux_g_nf(i1)
            ELSEIF(i1.lt.0) THEN
               i1=-i1
               flux_l=-flux_l_nf(i1)
               flux_g=-flux_g_nf(i1)
            ENDIF
         ELSE
            flux_l=0.d0
            flux_g=0.d0
         ENDIF
         IF(flux_l.gt.0.d0)THEN
            k=i
            flow_sp=flow_sp-cell%alphal(k)*cell%rhol(k)*flux_l
         ELSEIF(flux_l.lt.0.d0)THEN
            k=neigh(j+j0)
            flow_sp=flow_sp-cell%alphal(k)*cell%rhol(k)*flux_l
         ENDIF
         IF(flux_g.gt.0.d0)THEN
            k=i
            flow_sp=flow_sp-cell%alphag(k)*cell%rhog(k)*flux_g
         ELSEIF(flux_g.lt.0.d0)THEN
            k=neigh(j+j0)
            flow_sp=flow_sp-cell%alphag(k)*cell%rhog(k)*flux_g
         ENDIF
      ENDDO
      DO ix=1,num_f2
         i=nf2(1,ix)
         j0=i_neigh(i)-1
         j=nf2(2,ix)
         nf_number=get_nf_number(nbcon(j+j0))
         IF(nf_number.le.3) THEN
            CALL get_vector_disp(j,i,i1)
            IF(i1.gt.0) THEN
               flux_l=flux_l_nf(i1)
               flux_g=flux_g_nf(i1)
            ELSEIF(i1.lt.0) THEN
               i1=-i1
               flux_l=-flux_l_nf(i1)
               flux_g=-flux_g_nf(i1)
            ENDIF
         ELSE
            flux_l=0.d0
            flux_g=0.d0
         ENDIF
         IF(flux_l.gt.0.d0)THEN
            k=i
            flow_fd=flow_fd+cell%alphal(k)*cell%rhol(k)*flux_l
         ELSEIF(flux_l.lt.0.d0)THEN
            k=neigh(j+j0)
            flow_fd=flow_fd+cell%alphal(k)*cell%rhol(k)*flux_l
         ENDIF
         IF(flux_g.gt.0.d0)THEN
            k=i
            flow_fd=flow_fd+cell%alphag(k)*cell%rhog(k)*flux_g
         ELSEIF(flux_g.lt.0.d0)THEN
            k=neigh(j+j0)
            flow_fd=flow_fd+cell%alphag(k)*cell%rhog(k)*flux_g
         ENDIF
      ENDDO
      DO ix=1,num_f3
         i=nf3(1,ix)
         j0=i_neigh(i)-1
         j=nf3(2,ix)
         nf_number=get_nf_number(nbcon(j+j0))
         IF(nf_number.le.3) THEN
            CALL get_vector_disp(j,i,i1)
            IF(i1.gt.0) THEN
               flux_l=flux_l_nf(i1)
               flux_g=flux_g_nf(i1)
            ELSEIF(i1.lt.0) THEN
               i1=-i1
               flux_l=-flux_l_nf(i1)
               flux_g=-flux_g_nf(i1)
            ENDIF
         ELSE
            flux_l=0.d0
            flux_g=0.d0
         ENDIF
         IF(flux_l.gt.0.d0)THEN
            k=i
            flow_si=flow_si+cell%alphal(k)*cell%rhol(k)*flux_l
         ELSEIF(flux_l.lt.0.d0)THEN
            k=neigh(j+j0)
            flow_si=flow_si+cell%alphal(k)*cell%rhol(k)*flux_l
         ENDIF
         IF(flux_g.gt.0.d0)THEN
            k=i
            flow_si=flow_si+cell%alphag(k)*cell%rhog(k)*flux_g
         ELSEIF(flux_g.lt.0.d0)THEN
            k=neigh(j+j0)
            flow_si=flow_si+cell%alphag(k)*cell%rhog(k)*flux_g
         ENDIF
      ENDDO
!
      IF(np.gt.1)THEN
         tmp(1)=flow_sp
         tmp(2)=flow_fd
         CALL allreducei_r(tmp,2)
         flow_sp=tmp(1)
         flow_fd=tmp(2)
      ENDIF
!
      time1=time1+dt
      IF(time1.gt.0.1d0)THEN
         time1=0.d0
         volume=0.d0
         volume_sp=0.d0
         DO i=1,ncell_fluid
            IF(xloc(i,ndim).gt.0.d0) volume=volume+cell%alphal(i)*vol(i)
            IF(xloc(i,1).gt.1.04115.and.xloc(i,1).lt.1.38995.and.  &
               xloc(i,2).gt.1.04115.and.xloc(i,2).lt.1.38995)THEN
               IF(xloc(i,ndim).gt.0.25.and.xloc(i,ndim).lt.4.58)THEN
                  volume_sp=volume_sp+cell%alphal(i)*vol(i)
               ENDIF
            ENDIF
         ENDDO
!
         IF(np.gt.1)THEN
            tmp(1)=flow_si
            tmp(2)=volume
            tmp(3)=volume_sp
            CALL allreducei_r(tmp,3)
            flow_si  =tmp(1)
            volume   =tmp(2)
            volume_sp=tmp(3)
         ENDIF 
!
         IF(.false.)THEN  ! Plant test
            height=1.7546d0*(volume-8.2746d0)+9.549d0
         ELSE       ! VAPER test
            IF(volume.ge.30.093969d0)THEN
               area=(59.640839d0-30.093969d0)/5.d0
               height=5.545d0+(volume-30.093969d0)/area
            ELSEIF(volume.ge.7.132845d0)THEN
               area=(30.093969d0-7.132845d0)/4.d0
               height=1.545d0+(volume-7.132845d0)/area
            ELSEIF(volume.ge.3.746079d0)THEN
               area=(7.132845d0-3.746079d0)/0.59d0
               height=0.955d0+(volume-3.746079d0)/area
            ELSE
               area=3.746079d0/0.955d0
               height=volume/area       
            ENDIF
         ENDIF
!
         area_sp=0.12166144d0  ! 0.3488**2
         level_sp=volume_sp/area_sp
! 
         IF(n_top.gt.0)THEN
            p_top=p(n_top)
            t_top=cell%tg(n_top)
         ELSE
            p_top=-1.d10
            t_top=-1.d10
         ENDIF
!        IF(n_top_sp.gt.0)THEN
!           p_top_sp=p(n_top)
!        ELSE
!           p_top_sp=-1.d10
!        ENDIF
!        IF(n_top_fd.gt.0)THEN
!           p_top_fd=p(n_top)
!        ELSE
!           p_top_fd=-1.d10
!        ENDIF
!        IF(n_out1.gt.0)THEN
!           p_out1=p(n_top)
!        ELSE
!           p_out1=-1.d10
!        ENDIF
!        IF(n_out2.gt.0)THEN
!           p_out2=p(n_top)
!        ELSE
!           p_out2=-1.d10
!        ENDIF
         IF(np.gt.1)THEN
            tmp(1)=p_top
            tmp(2)=t_top
!           tmp(3)=p_top_sp
!           tmp(4)=p_top_fd
!           tmp(5)=p_out1
!           tmp(6)=p_out2
            CALL allreducei_max_r(tmp,2)
            p_top   =tmp(1) 
            t_top   =tmp(2) 
!           p_top_sp=tmp(3) 
!           p_top_fd=tmp(4) 
!           p_out1  =tmp(5) 
!           p_out2  =tmp(6) 
         ENDIF
!
         IF(myrank.eq.0)THEN
!           WRITE(999,90) time,flow_sp,flow_fd,flow_si,volume,height,level_sp,  &
!                          p_top,p_top_sp,p_top_fd,p_out1,p_out2,p_top-p_out2,   &
!                          p_top_fd-p_out1,t_top,k_sp,k_dp
            WRITE(333,90) time,flow_si,height,p_top/1000.d0,t_top
         ENDIF
      ENDIF
90    FORMAT(20(e16.9,1x))
!
      END SUBROUTINE fd_post_user
