!
      SUBROUTINE f_fric_face
!
      USE Zinterface
      USE Zparam       , ONLY: ndim
      USE Zcore        , ONLY: np
      USE Zvector      , ONLY: face_fr_l,face_fr_g
      USE Znum_cell    , ONLY: istart_nf
      USE Zvec_index   , ONLY: left_nf,right_non
      USE Zvec_geo     , ONLY: fac_non,fac1_non
      USE Zvector      , ONLY: vl_f_non,vg_f_non
      USE Zporous      , ONLY: l_mixing_vane
!
      IMPLICIT NONE
!     local variables
      INTEGER :: i
      INTEGER :: ii,kk
      INTEGER :: nf_number,istart,len,i1
!     local arrays
      REAL(8) :: vg_f_non1,vg_f_non2,vg_f_non3
      REAL(8) :: vl_f_non1,vl_f_non2,vl_f_non3
!
      IF(l_mixing_vane)then
         CALL udfn_f_fric_face
         RETURN
      ENDIF

      IF(np.gt.1) CALL communicate_2d(face_fr_l, &
                                      face_fr_g)
!
      IF(ndim.eq.2) THEN 
         nf_number=0
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
         DO i=1,len
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)
            vg_f_non1=face_fr_g(ii,1)*fac1_non(i)+face_fr_g(kk,1)*fac_non(i)
            vg_f_non2=face_fr_g(ii,2)*fac1_non(i)+face_fr_g(kk,2)*fac_non(i)
            vl_f_non1=face_fr_l(ii,1)*fac1_non(i)+face_fr_l(kk,1)*fac_non(i)
            vl_f_non2=face_fr_l(ii,2)*fac1_non(i)+face_fr_l(kk,2)*fac_non(i)
            vg_f_non(i1,1)=vg_f_non(i1,1)-vg_f_non1
            vg_f_non(i1,2)=vg_f_non(i1,2)-vg_f_non2
            vl_f_non(i1,1)=vg_f_non(i1,1)-vl_f_non1
            vl_f_non(i1,2)=vg_f_non(i1,2)-vl_f_non2
         ENDDO
!
      ELSE
!
         nf_number=0
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
         DO i=1,len
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)
            vg_f_non1=face_fr_g(ii,1)*fac1_non(i)+face_fr_g(kk,1)*fac_non(i)
            vg_f_non2=face_fr_g(ii,2)*fac1_non(i)+face_fr_g(kk,2)*fac_non(i)
            vg_f_non3=face_fr_g(ii,3)*fac1_non(i)+face_fr_g(kk,3)*fac_non(i)
            vl_f_non1=face_fr_l(ii,1)*fac1_non(i)+face_fr_l(kk,1)*fac_non(i)
            vl_f_non2=face_fr_l(ii,2)*fac1_non(i)+face_fr_l(kk,2)*fac_non(i)
            vl_f_non3=face_fr_l(ii,3)*fac1_non(i)+face_fr_l(kk,3)*fac_non(i)
            vg_f_non(i1,1)=vg_f_non(i1,1)-vg_f_non1
            vg_f_non(i1,2)=vg_f_non(i1,2)-vg_f_non2
            vg_f_non(i1,3)=vg_f_non(i1,3)-vg_f_non3
            vl_f_non(i1,1)=vl_f_non(i1,1)-vl_f_non1
            vl_f_non(i1,2)=vl_f_non(i1,2)-vl_f_non2
            vl_f_non(i1,3)=vl_f_non(i1,3)-vl_f_non3
         ENDDO
!
      ENDIF
!
      RETURN
      END SUBROUTINE f_fric_face
!
!======================================================================
!======================================================================
!
      SUBROUTINE udfn_f_fric_face
!
      USE Zinterface
      USE VOL_DATA     , ONLY: cell
      USE Zcore        , ONLY: np
      USE Znum_cell    , ONLY: istart_nf,i_neigh
      USE Zvec_index   , ONLY: left_nf,right_non,jneigh_nf
      USE Zvec_param   , ONLY: nf_non
      USE Zporous      , ONLY: mv_fac,mv_non
      USE Zvector      , ONLY: vl_f_non
      USE Zvec_geo     , ONLY: fac_non,fac1_non

      IMPLICIT NONE

      INTEGER :: i
      INTEGER :: ii,kk
      INTEGER :: nf_number,istart,len,i1
      INTEGER :: j,j0
      
      LOGICAL,SAVE :: initial
      DATA initial /.true./
!
      IF(np.gt.1) CALL communicate_1d(cell%vfwl_x, &
                                      cell%vfwl_y, &
                                      cell%vfwl_z)

         IF(initial)then
            initial=.false.

            ALLOCATE(mv_non(nf_non))
            mv_non=0
            
            nf_number=0
            istart=istart_nf(1,nf_number)
            len   =istart_nf(2,nf_number)
            DO i=1,len
               i1=istart+i
               ii=left_nf(i1)
               j0=i_neigh(ii)-1
               j=jneigh_nf(i1)+j0
               mv_non(i1) =mv_fac(j)
            ENDDO
         ENDIF

         nf_number=0
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
         DO i=1,len
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)     
            vl_f_non(i1,1)=cell%vfwl_x(ii)*fac1_non(i)*mv_non(i1)+ &
                           cell%vfwl_x(kk)*fac_non(i) *mv_non(i1)
            vl_f_non(i1,2)=cell%vfwl_y(ii)*fac1_non(i)*mv_non(i1)+ &
                           cell%vfwl_y(kk)*fac_non(i) *mv_non(i1)
            vl_f_non(i1,3)=dmin1(cell%vfwl_z(ii),cell%vfwl_z(kk))
         ENDDO

      RETURN
      END SUBROUTINE udfn_f_fric_face
