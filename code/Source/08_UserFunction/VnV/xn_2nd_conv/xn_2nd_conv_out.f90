!
      SUBROUTINE xn_2nd_conv_out
!
      USE VOL_DATA     , ONLY: cell
      USE Zzone        , ONLY: ncell_fluid,ncell_fluid_all
      USE Zcore        , ONLY: myrank
      USE Ztimecon     , ONLY: time
      USE Zcoord1      , ONLY: xloc
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i,ii,na
      LOGICAL,SAVE :: initial=.true.
!.....Local arrays
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: dat
!
      IF(initial) THEN
         IF(myrank.eq.0) THEN
            OPEN(990,file='VFT14_xn_2nd_conv_ref.dat')
         ENDIF
         initial=.false.
      ENDIF
!
      IF(time.ge.1000.d0)THEN
         na=ncell_fluid_all
         IF(myrank.eq.0) THEN
            ALLOCATE(dat(na,2))
         ELSE
            ALLOCATE(dat(1,2))
         ENDIF
         CALL gatherv_r(xloc(1,3) ,ncell_fluid,dat(1,1),na,0)
         CALL gatherv_r(cell%quala,ncell_fluid,dat(1,2),na,0)
!
         IF(myrank.eq.0) THEN
            DO i=1,48
               ii=11*(i-1)+6
               WRITE(990,10) dat(ii,1),dat(ii,2)
            ENDDO
         ENDIF
      ENDIF
!
   10 FORMAT(2(4e15.7,1x))
!
      END SUBROUTINE xn_2nd_conv_out
