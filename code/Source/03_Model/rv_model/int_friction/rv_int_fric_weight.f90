      REAL(8) FUNCTION rv_int_fric_weight(alphag, alph_left, alph_right)
!
!.....Weighting factor by an independent variable's position and its interval.
!     When alphag is within its left and right boundary, the weighting factor is zero.
!     If alphag goes beyond its interval, the weighting factor linearly fades away to
!     zero within the fading width vt_l and vt_r.
!
      IMPLICIT NONE
!      
      REAL(8),INTENT(IN) :: alphag       ! a variable.
      REAL(8),INTENT(IN) :: alph_left    ! left limit for the variable.
      REAL(8),INTENT(IN) :: alph_right   ! right limit for the variable.
!     REAL(8),INTENT(IN) :: vtran_max    ! transition width
!
!     REAL(8),PARAMETER :: vtran_min = 0.01 ! transition width
!     REAL(8),PARAMETER :: vtran_max = 0.01 ! transition width
      REAL(8),PARAMETER :: vtran_width=100.d0   !< transition width
!
!     REAL(8)::wl,wr,vt_l,vt_r
      REAL(8)::wl,wr
   
!     IF(alph_left<=0.0.or.1.0<=alph_left)THEN
!        vt_l=vtran_min
!     ELSE
!        vt_l=vtran_max
!     ENDIF
!  
!     IF(alph_right<=0.0.or.1.0<=alph_right)THEN
!        vt_r=vtran_min
!     ELSE
!        vt_r=vtran_max
!     ENDIF
   
!     wl=0.5+0.5*(alphag-alph_left)/vt_l
      wl=0.5+0.5*(alphag-alph_left)*vtran_width
      wl=DMAX1(0.0,DMIN1(wl,1.0))
!     wr=0.5+0.5*(alphag-alph_right)/vt_r
      wr=0.5+0.5*(alphag-alph_right)*vtran_width
      wr=DMAX1(0.0,DMIN1(wr,1.0))
      rv_int_fric_weight=wl-wr
!      
      RETURN
      ENDFUNCTION rv_int_fric_weight
!--------------------------------------------------------------------------------------
      SUBROUTINE rv_int_fric_weight_interval(x,x1,x2,wf1,wf2)
!
!.....Check if wf1 and wf2 are continuous. They jump at x1 and x2.
!
      IMPLICIT NONE
!      
      REAL(8) x,x1,x2,wf1,wf2,dx
      IF(x.le.x1)THEN
         wf1=1.0d0
         wf2=0.0d0
      ELSEIF(x.ge.x2)THEN
         wf1=0.0d0
         wf2=1.0d0
      ELSE
         dx=x2-x1
         wf1=(x-x1)/dx
         wf2=(x2-x)/dx
      ENDIF
!
      RETURN
      ENDSUBROUTINE rv_int_fric_weight_interval
