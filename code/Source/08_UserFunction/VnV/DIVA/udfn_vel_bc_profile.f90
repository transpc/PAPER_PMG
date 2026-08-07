!
      SUBROUTINE udfn_vel_bc_profile(icell,nin,f_profile)
!
!     This routine defines the velocity disribution of inlet boundary condition
!
      USE Zconst1         , ONLY: vv_prob       
      USE Zcoord1         , ONLY: xloc
      USE Zmpi            , ONLY: jperm      
      USE Zio_unit        , ONLY: unit_log
!
      IMPLICIT NONE
!
      INTEGER icell,nin,i,j,profile_direc
      INTEGER,SAVE:: index_cellnum(16,16)    !ce15x15 (LSJ)
      REAL(8) x,f_profile
      REAL(8) a1,a2,a3,a4,loc_min,loc_max,length_tot
      REAL(8),SAVE:: f_profile_tmp(16,16)    !ce15x15 (LSJ)
!
      DATA profile_direc /3/                                                ! select x,y,z dimension 
      DATA a1,a2,a3,a4 /0.63401d0, 2.78251d0,-5.67148d0,3.14998d0/          ! coeff. of poly line (4th)
      DATA loc_min,loc_max /0.177419d0,0.822581d0/                          ! min, max location to be applied (normalized)
      DATA length_tot /1.0d0/                                              ! totla length to normalize  (org=0.62d0)
!
      LOGICAL,SAVE::INITIAL
      DATA INITIAL/.TRUE./
!  
      IF(vv_prob.eq.'ce15x15')THEN  
         IF(initial) THEN
            initial=.false.
            OPEN(111,file='index_ce_inlet.dat')
            OPEN(112,file='index_ce_f_profile.dat')
            DO i=1,16
               READ(111,*) (index_cellnum(i,j),j=1,16)
               READ(112,*) (f_profile_tmp(i,j),j=1,16)
            ENDDO
         ENDIF
!         
         DO i=1,16
         DO j=1,16
            IF(jperm(icell).eq.index_cellnum(i,j))THEN
               f_profile=f_profile_tmp(i,j)
            ENDIF   
         ENDDO
         ENDDO   
      
      ELSEIF(vv_prob.eq.'DIVA-NEW')THEN
         IF(nin.eq.1)THEN                                                      ! Select velocity B.C to be applied (1 to 4)
            x=xloc(icell,profile_direc)/length_tot
            IF(x.le.loc_min)THEN
               x=loc_min
               f_profile=a1+a2*x+a3*x**2+a4*x**3
            ELSEIF(x.gt.loc_max)THEN
               x=loc_max
               f_profile=a1+a2*x+a3*x**2+a4*x**3
            ELSE
               f_profile=a1+a2*x+a3*x**2+a4*x**3
            ENDIF
         ELSE
            f_profile=1.0d0
         ENDIF    
!             
         IF(f_profile.le.0.0d0.or.f_profile.ge.2.0d0)THEN
           WRITE(*,*) 'fprofile=',f_profile,'x=',x
           WRITE(unit_log,*) 'fprofile=',f_profile,'x=',x
           STOP
         ENDIF     
      ELSE
         f_profile=1.0d0
      ENDIF   
!      
      END SUBROUTINE udfn_vel_bc_profile 
!
      SUBROUTINE udfn_vel_bc_profile_inl
!
!     This routine defines the velocity disribution of inlet boundary condition
!
      USE Zmpi         , ONLY: jperm      
      USE Zcore        , ONLY: np,myrank
      USE Znum_cell    , ONLY: istart_nf,istart_nbcon_nf
      USE Zvec_index   , ONLY: left_nf,nbcon_nf
      USE Zcoord1      , ONLY: xloc
      USE Zconst1      , ONLY: vv_prob       
      USE Zuserdefined , ONLY: vel_bc_profile_inl
      USE Zio_unit     , ONLY: unit_log
!
      IMPLICIT NONE
!
      INTEGER :: icell,i,j,profile_direc
      INTEGER :: ii,jj,k
      INTEGER :: nf_number,istart,len,istart2,i1,i2
      INTEGER :: kill
      INTEGER :: index_cellnum(16,16)    !ce15x15 (LSJ)
      REAL(8) :: x,f_profile
      REAL(8) :: a1,a2,a3,a4,loc_min,loc_max,length_tot
      REAL(8) :: f_profile_tmp(16,16)    !ce15x15 (LSJ)
!
      DATA profile_direc /3/                                                ! select x,y,z dimension 
      DATA a1,a2,a3,a4 /0.63401d0, 2.78251d0,-5.67148d0,3.14998d0/          ! coeff. of poly line (4th)
      DATA loc_min,loc_max /0.177419d0,0.822581d0/                          ! min, max location to be applied (normalized)
      DATA length_tot /1.0d0/                                               ! totla length to normalize  (org=0.62d0)
!
      IF(vv_prob.eq.'ce15x15')THEN  
         IF(myrank.eq.0) THEN
            OPEN(111,file='index_ce_inlet.dat')
            OPEN(112,file='index_ce_f_profile.dat')
            DO i=1,16
               READ(111,*) (index_cellnum(i,j),j=1,16)
               READ(112,*) (f_profile_tmp(i,j),j=1,16)
            ENDDO
            CLOSE(111)
            CLOSE(112)
         ENDIF
         IF(np.gt.1) THEN
            CALL broadcast_i(index_cellnum,16*16)
            CALL broadcast_r(f_profile_tmp,16*16)
         ENDIF
!         
         nf_number=2
         istart=istart_nf(1,nf_number)
         istart2=istart_nbcon_nf(nf_number)
         len    =istart_nf(2,nf_number)
         DO i=1,len
            i1=istart+i
            i2=istart2+i
            icell=left_nf(i1)
            k=nbcon_nf(i2)
            DO jj=1,16
               DO ii=1,16
                  IF(jperm(icell).eq.index_cellnum(ii,jj))THEN
                     f_profile=f_profile_tmp(ii,jj)
                  ENDIF   
               ENDDO
            ENDDO   
            vel_bc_profile_inl(i)=f_profile
         ENDDO   
      ELSEIF(vv_prob.eq.'DIVA-NEW')THEN
         nf_number=2
         istart=istart_nf(1,nf_number)
         istart2=istart_nbcon_nf(nf_number)
         len    =istart_nf(2,nf_number)
         kill=0
         DO i=1,len
            i1=istart+i
            i2=istart2+i
            icell=left_nf(i1)
            k=nbcon_nf(i2)
!
            IF(k.eq.1)THEN                                                      ! Select velocity B.C to be applied (1 to 4)
               x=xloc(icell,profile_direc)/length_tot
               IF(x.le.loc_min)THEN
                  x=loc_min
                  f_profile=a1+a2*x+a3*x**2+a4*x**3
               ELSEIF(x.gt.loc_max)THEN
                  x=loc_max
                  f_profile=a1+a2*x+a3*x**2+a4*x**3
               ELSE
                  f_profile=a1+a2*x+a3*x**2+a4*x**3
               ENDIF
               IF(f_profile.le.0.0d0.or.f_profile.ge.2.0d0)THEN
                  WRITE(*,*) myrank,'fprofile=',f_profile,'x=',x
                  WRITE(unit_log,*) 'fprofile=',f_profile,'x=',x
                  kill=1
                  exit
               ENDIF     
            ELSE
               f_profile=1.0d0
            ENDIF    
            vel_bc_profile_inl(i)=f_profile
         ENDDO
         IF(np.gt.1) CALL allreducei_i1(kill)
         IF(kill.gt.0) THEN
            CALL barrier_mpi
            CALL finalize_mpi
            STOP
         ENDIF
      ELSE
         nf_number=2
         len    =istart_nf(2,nf_number)
         DO i=1,len
            vel_bc_profile_inl(i)=1.0d0
         ENDDO
      ENDIF   
!      
      END SUBROUTINE udfn_vel_bc_profile_inl
      
