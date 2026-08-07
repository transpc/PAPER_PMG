!
      SUBROUTINE udfn_subchannel_geo(nn,vol,poro,hd,sl,sgap)
!
!     Define user-defined geometry data for subchannel-scale applications
!
      USE Zcore       , ONLY: myrank,np
      USE Zconst1     , ONLY: vv_prob
      USE Zconst2     , ONLY: hydraulicd_init
      USE Zporous     , ONLY: chn_type_tmp
      USE Zrv_ncell   , ONLY: asm_nx,asm_ny,asm_nz,asm_ni,asm_ni2,chn_nx,chn_ny,chn_nz
      USE Zparam      , ONLY: ndim
      USE Znum_cell   , ONLY: i_neigh_tmp,sv_tmp1,perm_tmp1

      USE Zporous     , ONLY: n_sg,h_sg,sg_loc
      USE Zporous     , ONLY: l_spacer_grid
      USE Zio_unit    , ONLY: unit_log
!
      IMPLICIT NONE

      !input
      INTEGER nn
      REAL(8) vol(nn)
      !output
      REAL(8) poro(nn),hd(nn)
      REAL(8) perm(nn,ndim),sl(nn,ndim),sgap(nn,ndim)

      !open file index
      INTEGER unit_sc,scin      
      
      INTEGER i,ii,j
      INTEGER i0,z0,iasm
      INTEGER i1,i2,i3,i4,i5,i6,i7
      REAL(8) sa
      REAL(8) r1,r2,r3,r4,r5,r6,r7,r8
!     permeability
      REAL(8) xn_x,xn_y,xn_z


       DATA scin/1/
!   
!======================================================================
! READ Subchannel geometric information
!======================================================================

      unit_sc=15501
      OPEN(unit=unit_sc,file='subchannel_info.in',status='old', iostat=scin)
      IF(scin.ne.0)then
         IF(myrank.eq.0) WRITE(*       ,*)'Program was terminated due to lack of <subchannel_info.in>!'
         IF(myrank.eq.0) WRITE(unit_log,*)'Program was terminated due to lack of <subchannel_info.in>!'
         CALL finalize_mpi
      ENDIF
      
      READ(unit_sc,*)i0,z0,iasm
      i0=i0*z0

      ALLOCATE(chn_nx(nn),chn_ny(nn),chn_nz(nn))
      ALLOCATE(chn_type_tmp(nn))
      
      DO i=1,nn
         chn_nx(i)      =0
         chn_ny(i)      =0
         chn_nz(i)      =0          !PSH
         chn_type_tmp(i)=0
      ENDDO
      IF(myrank.eq.0)then
      DO i=1,nn
         poro(i)        =1.0d0
         perm(i,:)      =1.0d0
         hd(i)          =1.0d0
         sl(i,:)        =1.0d0
         sgap(i,:)      =0.0d0
      ENDDO
      ENDIF
      
!.....For Multi-Assembly application      
      IF(iasm.eq.1)then
         ALLOCATE(asm_nx(nn),asm_ny(nn),asm_nz(nn),asm_ni(nn),asm_ni2(nn))
         DO i=1,nn
            asm_nx(i)=0
            asm_ny(i)=0
            asm_nz(i)=0
            asm_ni(i)=0
            asm_ni2(i)=0
         ENDDO
      ENDIF   

!.....Read subchannel information      

!.....Single assembly application          
      IF(myrank.eq.0)then
      IF(iasm.eq.0)then 
         DO ii=1,i0
            READ(unit_sc,*)i,i1,i2,i3,i4,r1,r2,r3,r4,r5,r6,r7,r8
            chn_nx(i)=i1
            chn_ny(i)=i2
            chn_nz(i)=i3
            chn_type_tmp(i)=i4
            poro(i)  =r1
            perm(i,1)=r2
            perm(i,2)=r3
            perm(i,3)=MAX(0.3,r1)
            hd(i)    =r4
            sl(i,1)  =r5
            sl(i,2)  =r6
            sgap(i,1)=r7
            sgap(i,2)=r8
         ENDDO
!
!.....Multi-assembly (TO BE modified later)            
      ELSE
         DO ii=1,i0
            READ(unit_sc,*)i,i1,i2,i3,i4,i5,i6,i7,r1,r2,r3,r4,r5,r6,r7,r8
            asm_nx(i)=i1
            asm_ny(i)=i2
            asm_nz(i)=i3
            chn_nx(i)=i4
            chn_ny(i)=i5
            chn_nz(i)=i6
            chn_type_tmp(i)=i7
            poro(i)  =r1
            perm(i,1)=r2
            perm(i,2)=r3
            perm(i,3)=MAX(0.3,r1)
            hd(i)    =r4
            sl(i,1)  =r5
            sl(i,2)  =r6
            sgap(i,1)=r7
            sgap(i,2)=r8
         ENDDO
      ENDIF
!
!.....Permeability setting
!
      IF(ndim.eq.2) THEN
         DO i=1,nn
            DO j=i_neigh_tmp(i),i_neigh_tmp(i+1)-1
               sa=sv_tmp1(j,1)**2+sv_tmp1(j,2)**2
               sa=1.d0/DSQRT(sa)
               xn_x=sv_tmp1(j,1)*sa 
               xn_y=sv_tmp1(j,2)*sa
               IF(ABS(xn_x).gt.0.1) perm_tmp1(j)=perm(i,1)
               IF(ABS(xn_y).gt.0.1) perm_tmp1(j)=perm(i,2)
            ENDDO
         ENDDO
       ELSE
         DO i=1,nn
            DO j=i_neigh_tmp(i),i_neigh_tmp(i+1)-1
               sa=sv_tmp1(j,1)**2+sv_tmp1(j,2)**2+sv_tmp1(j,3)**2
               sa=1.d0/DSQRT(sa)
               xn_x=sv_tmp1(j,1)*sa
               xn_y=sv_tmp1(j,2)*sa
               xn_z=sv_tmp1(j,3)*sa
               IF(ABS(xn_x).gt.0.1) perm_tmp1(j)=perm(i,1)
               IF(ABS(xn_y).gt.0.1) perm_tmp1(j)=perm(i,2)
               IF(ABS(xn_z).gt.0.1) perm_tmp1(j)=perm(i,3)
            ENDDO
         ENDDO
      ENDIF
      
      ENDIF !IF(myrank.~

      hydraulicd_init=1   
      
      IF(np.gt.1) THEN
         CALL broadcast_i(chn_nx,nn)
         CALL broadcast_i(chn_ny,nn)
         CALL broadcast_i(chn_nz,nn)
         CALL broadcast_i(chn_type_tmp,nn)
         IF(iasm.eq.1)then      !PSH
            CALL broadcast_i(asm_nx,nn)
            CALL broadcast_i(asm_ny,nn)
            CALL broadcast_i(asm_nz,nn)
            CALL broadcast_i(asm_ni,nn)
            CALL broadcast_i(asm_ni2,nn)
         ENDIF    
      ENDIF
      
!
!........Spacer Grid
!
      IF(l_spacer_grid)then
         IF(myrank.eq.0)then 
            READ(unit_sc,*) n_sg
            IF(n_sg.eq.0)then
               WRITE(unit_log,*)'          Spacer Grid is not defined in subchannel_info.in'
               WRITE(*       ,*)'          Spacer Grid is not defined in subchannel_info.in'
               goto 1551
            ENDIF
            ALLOCATE(h_sg(n_sg),sg_loc(n_sg))
            DO i=1,n_sg
               READ(unit_sc,*)h_sg(i),sg_loc(i)
            ENDDO
         ENDIF   
         IF(np.gt.1) THEN
            CALL broadcast_i1(n_sg)
         ENDIF
         IF(myrank.ne.0) then
             ALLOCATE(  h_sg(n_sg))
             ALLOCATE(sg_loc(n_sg))
         ENDIF    
         IF(np.gt.1) CALL broadcast_r(  h_sg,n_sg)
         IF(np.gt.1) CALL broadcast_r(sg_loc,n_sg)

         IF(myrank.eq.0) then
            WRITE(unit_log,*)'          Spacer Grid is successfully read from subchannel_info.in'
            WRITE(*       ,*)'          Spacer Grid is successfully read from subchannel_info.in'
         ENDIF
1551  continue
      ELSE
         IF(myrank.eq.0)then
            WRITE(unit_log,*)'          Spacer Grid is not defined in subchannel_info.in'
            WRITE(*       ,*)'          Spacer Grid is not defined in subchannel_info.in'
         ENDIF
      ENDIF
!
      CLOSE(unit_sc)

!
!.....Specific geometric modification for each examples
!
      IF(myrank.eq.0)then
          
         IF(vv_prob.eq.'GE3x3') then     !!!CYJ_tmp       
            DO i=1,nn
               IF(chn_type_tmp(i).eq.3)then
                  vol(i)=vol(i)-0.000001220808265 
                  DO j=i_neigh_tmp(i),i_neigh_tmp(i+1)-1
                     IF(sv_tmp1(j,3).gt.0.d0)THEN
                        sv_tmp1(j,3)=sv_tmp1(j,3)-0.000024518d0
                     ELSEIF(sv_tmp1(j,3).lt.0.d0)THEN
                        sv_tmp1(j,3)=sv_tmp1(j,3)+0.000024518d0
                     ENDIF
                  ENDDO   
               ENDIF
            ENDDO
         ENDIF   
         
      ENDIF

      RETURN
      END SUBROUTINE udfn_subchannel_geo
!
