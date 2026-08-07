!
      SUBROUTINE multi_ncg_out
!
      USE VOL_DATA     , ONLY: cell
      USE Zcore        , ONLY: myrank
      USE Zncg         , ONLY: qn_cell 
      USE Ztimecon     , ONLY: time
      USE Zzone        , ONLY: ncell_fluid,ncell_fluid_all
!
      IMPLICIT NONE
!
      INTEGER,PARAMETER :: nout=10
!.....Local variables
      INTEGER :: na
      INTEGER,SAVE :: iout=0
      LOGICAL,SAVE :: initial=.true.
!.....Local arrays
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: dat
!
      IF(initial.and.myrank.eq.0)THEN
         OPEN(990,file='VFT12_multi_ncg_ref.dat')
         initial=.false.
      ENDIF
!
      iout=iout+1
      IF(iout.ge.nout)THEN
         na=ncell_fluid_all
         IF(myrank.eq.0) THEN
            ALLOCATE(dat(na,3))
         ELSE
            ALLOCATE(dat(1,3))
         ENDIF
         CALL gatherv_r(cell%tg     ,ncell_fluid,dat(1,1),na,0)
         CALL gatherv_r(qn_cell(1,1),ncell_fluid,dat(1,2),na,0)
         CALL gatherv_r(qn_cell(1,2),ncell_fluid,dat(1,3),na,0)
         IF(myrank.eq.0) THEN
            WRITE(990,10) time,dat(55,1),dat(955,1),dat(55,2),dat(955,2),dat(55,3),dat(955,3)
         ENDIF
         iout=0
         DEALLOCATE(dat)
      ENDIF
!
   10 FORMAT(10(4e15.7,1x))
!
      END SUBROUTINE multi_ncg_out
