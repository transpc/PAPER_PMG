!
      SUBROUTINE grad_press_ls_c_2d
!
!     This routine calculate geometry related matrix coefficient
!     for 2D least square method.
!
      USE Zzone        , ONLY: ncell_fluid
      USE Zparam       , ONLY: ndim      
      USE Zvec_param   , ONLY: nf_non
      USE Znum_cell    , ONLY: istart_nf
      USE Zvec_index   , ONLY: left_nf,right_non
      USE Zcoord1      , ONLY: xloc
      USE Zgrad_ls_c2d , ONLY: a_2,a11_2,a12_2,a22_2,det_2
      USE Zbc_index    , ONLY: ngrad
      USE Zgrad_ls_c3d , ONLY: lsindex
      USE Zrv_model    , ONLY: rv_valve     
!
      IMPLICIT NONE
!.....Local variables
      INTEGER :: i
      INTEGER :: ii,kk
      INTEGER :: nf_number,len,istart,i1
      REAL(8) :: dx1,dx2
      REAL(8) :: b11
      REAL(8) :: b12,b22
      REAL(8) :: eps
!.....Local vector arrays
      REAL(8) :: b_non1(nf_non,ndim)
      REAL(8) :: b_non2(nf_non,ndim:ndim)
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
         IF(     (ngrad(ii).gt.0 .and. lsindex(ii).ne.0)       &
            .or. (ngrad(kk).gt.0 .and. lsindex(kk).ne.0)) THEN
            dx1=xloc(kk,1)-xloc(ii,1)
            dx2=xloc(kk,2)-xloc(ii,2)
            b_non1(i,1)=dx1*dx1
            b_non1(i,2)=dx1*dx2
            b_non2(i,2)=dx2*dx2
         ENDIF
      ENDDO
!      
!.....valve model         
!     
      IF(rv_valve.eq.1) CALL valve_model_grad_press_ls_c_2d(b_non1,b_non2)
!
      CALL grad_press_ls_c_2d_sum(b_non1,b_non2, &
                                  a11_2,a12_2,a22_2)
!
      DO i=1,ncell_fluid
         IF(ngrad(i).gt.0 .and. lsindex(i).ne.0) THEN
            b11=a11_2(i)
            b12=a12_2(i)
            b22=a22_2(i)
            IF(b11.lt.eps) b11=1.0d0
            IF(b22.lt.eps) b22=1.0d0
            a11_2(i)=b22
            a12_2(i)=-b12
            a22_2(i)=b11
            a_2(i)=b11*b22-b12*b12
            IF(abs(a_2(i)).gt.eps) det_2(i)=1.0d0
         ENDIF
      ENDDO
!
      END SUBROUTINE grad_press_ls_c_2d
