!
      SUBROUTINE grad_press_ls_c_3d
!
!     This routine calculate geometry related matrix coefficient
!     for 3D least square method.
!
      USE Zzone        , ONLY: ncell_fluid
      USE Zparam       , ONLY: ndim      
      USE Zvec_param   , ONLY: nf_non
      USE Znum_cell    , ONLY: istart_nf
      USE Zvec_index   , ONLY: left_nf,right_non
      USE Zcoord1      , ONLY: xloc
      USE Zgrad_ls_c3d , ONLY: a_3,a11_3,a12_3,a22_3,a13_3,a23_3,a33_3,det_3,lsindex
      USE Zbc_index    , ONLY: ngrad
      USE Zrv_model    , ONLY: rv_valve         
!
      IMPLICIT NONE
!
!     local variables
      INTEGER :: i
      INTEGER :: ii,kk
      INTEGER :: nf_number,len,istart,i1
      REAL(8) :: dx1,dx2,dx3
      REAL(8) :: b11
      REAL(8) :: b12,b22
      REAL(8) :: b13,b23,b33
      REAL(8) :: s1,s2,s3
      REAL(8) eps
!     local arrays
      REAL(8) :: b_non1(nf_non,ndim)
      REAL(8) :: b_non2(nf_non,2:ndim)
      REAL(8) :: b_non3(nf_non,ndim:ndim)
!
      DATA eps/1.e-15/
! 
!.....Cells non
!
      nf_number=0
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      DO i=1,len  
         i1=istart+i
         ii=left_nf(i1)
         kk=right_non(i)
         IF(     (ngrad(ii).gt.0 .and. lsindex(ii).ne.0)            &
            .or. (ngrad(kk).gt.0 .and. lsindex(kk).ne.0)) THEN
            dx1=xloc(kk,1)-xloc(ii,1)
            dx2=xloc(kk,2)-xloc(ii,2)
            dx3=xloc(kk,3)-xloc(ii,3)
            b_non1(i,1)=dx1*dx1
            b_non1(i,2)=dx1*dx2
            b_non1(i,3)=dx1*dx3
            b_non2(i,2)=dx2*dx2
            b_non2(i,3)=dx2*dx3
            b_non3(i,3)=dx3*dx3
         ENDIF
      ENDDO
!         
!.....valve model         
!       
      IF(rv_valve.eq.1) CALL valve_model_grad_press_ls_c_3d(b_non1,b_non2,b_non3)
!
      CALL grad_press_ls_c_3d_sum(b_non1,b_non2,b_non3, &
                                  a11_3,a12_3,a13_3,a22_3,a23_3,a33_3)
!
      DO i=1,ncell_fluid
         IF(ngrad(i).gt.0 .and. lsindex(i).ne.0) THEN
            b11=a11_3(i)
            b12=a12_3(i)
            b13=a13_3(i)
            b22=a22_3(i)
            b23=a23_3(i)
            b33=a33_3(i)
            IF(b11.lt.eps) b11=1.0d0
            IF(b22.lt.eps) b22=1.0d0
            IF(b33.lt.eps) b33=1.0d0
            s1=b22*b33-b23*b23 !L11
            s2=b13*b23-b12*b33 !L12
            s3=b12*b23-b13*b22 !L13
            a11_3(i)=s1
            a12_3(i)=s2
            a13_3(i)=s3
            a22_3(i)=b11*b33-b13*b13     !L22
            a23_3(i)=b13*b12-b11*b23     !L23
            a33_3(i)=b11*b22-b12*b12     !L33
            a_3(i)= b11*s1+b12*s2+b13*s3 !L
            IF(abs(a_3(i)).gt.eps) det_3(i)=1.0d0  ! meaningful a(i)=L
         ENDIF
      ENDDO
!
      END SUBROUTINE grad_press_ls_c_3d
