MODULE unitManager
   IMPLICIT NONE
   INTEGER :: cnt = 0, next
   INTEGER, PARAMETER :: offset = 1000
   INTEGER, PARAMETER :: shift = 10000    !< Translate starting index for avoiding conflict with other codes.
   INTEGER :: units(offset)
   CHARACTER(LEN=40) :: description(offset)
 !  
CONTAINS
   FUNCTION createUnit(dsc)
      USE Zcore, ONLY: myrank
      USE Zio_unit, ONLY: unit_keyborad, unit_screen
      IMPLICIT NONE
      INTEGER :: createUnit, i
      LOGICAL :: check
!     CHARACTER(20) :: name
      CHARACTER(*), OPTIONAL :: dsc
      next=cnt+1
      DO i = next+offset*myrank, offset+offset*myrank
         createUnit = i
         INQUIRE(unit=createUnit,opened=check)
         IF(check==.false.) THEN
            cnt = i - offset*myrank
            EXIT
         ENDIF
      ENDDO

      IF ( cnt > offset ) THEN
         PRINT *, "ERROR :: Too many units are USEd in myrank : ", myrank
         createUnit = 0
         RETURN
      END IF
      
      createUnit = createUnit+shift

      units(cnt)=createUnit
      IF(present(dsc)) THEN
         description(cnt) = dsc
      ELSE
         description(cnt) = 'no description'  
      ENDIF
   END FUNCTION createUnit

   FUNCTION getUnit(dsc)
      USE Zcore, ONLY: myrank
      IMPLICIT NONE
      CHARACTER(*) :: dsc
      INTEGER :: i, getUnit
      LOGICAL :: find = .false.
      DO i = 1, offset
         IF ( trim(description(i)) == trim(dsc) ) THEN
            find = .true.
            EXIT
         ENDIF
      END DO
      IF(find) THEN
         getUnit=myrank*offset+i
      ELSE
         PRINT *, "Unit Error :: '", dsc, "' is not found"
         getUnit=0
      ENDIF
   END function getUnit
!   
!...For debugging.
   subroutine PRINTAllUnits
      USE Zcore, ONLY: myrank, np
      USE Zio_unit, ONLY: unit_keyborad, unit_screen
      implicit none
!DEC$IF defined (mpi_flag)
      INCLUDE 'mpif.h'
!DEC$ENDIF      
      INTEGER :: i, ierr
      INTEGER :: cnt_all, cnt_arr(np), cnt_dsp(np)
      CHARACTER(LEN=40), ALLOCATABLE :: description_all(:)
      INTEGER, ALLOCATABLE :: units_all(:)

      IF(myrank == 0) then
         units(unit_keyborad) = unit_keyborad
         units(unit_screen) = unit_screen
         description(unit_keyborad)="RESERVED for keyboard"
         description(unit_screen)="RESERVED for screen"
      ENDIF
!DEC$IF defined (mpi_flag)
      cnt_arr=0
      CALL MPI_GATHER(cnt    ,1,MPI_INTEGER,&
                     cnt_arr,1,MPI_INTEGER,&
                     0,MPI_COMM_WORLD,ierr)
      cnt_dsp = 0
      DO i = 2, np
         cnt_dsp(i) = cnt_dsp(i-1)+cnt_arr(i-1)
      END DO
      cnt_all = sum(cnt_arr)
      ALLOCATE(description_all(cnt_all), units_all(cnt_all))
      CALL MPI_GATHERV(description    ,cnt*40               ,MPI_CHARACTER,&
                     description_all,cnt_arr*40,cnt_dsp*40,MPI_CHARACTER,&
                     0,MPI_COMM_WORLD,ierr)
      CALL MPI_GATHERV(units    ,cnt            ,MPI_INTEGER,&
                     units_all,cnt_arr,cnt_dsp,MPI_INTEGER,&
                     0,MPI_COMM_WORLD,ierr)
      IF ( myrank == 0 ) then
         DO i = 1, cnt_all
            PRINT *, units_all(i), description_all(i) 
         END DO   
      END IF
!DEC$ELSE
      DO i = 1, cnt
         PRINT *, units(i), description(i) 
      END DO   
!DEC$ENDIF
       
   END SUBROUTINE PRINTAllUnits
END MODULE unitManager
