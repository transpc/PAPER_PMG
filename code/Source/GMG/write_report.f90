     SUBROUTINE STOP_MPI(report_text)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
!    STOP MPI and writing a report to fort.999
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
      USE MD_MPI, ONLY: myrank
	  
      IMPLICIT NONE
	  
!DEC$IF defined (mpi_flag)
      include 'mpif.h'
!DEC$ENDIF

      CHARACTER :: report_text*100
	  
      INTEGER(4) rc
	  IF(myrank.EQ.0) THEN
	    WRITE(999,*) report_text
	    WRITE(*,*) 'PMG STOP: ', report_text
	  ENDIF
	  
! 
!DEC$IF defined (mpi_flag)
      call mpi_finalize(rc)
!DEC$ENDIF

      STOP
	  
	  RETURN
	  
      ENDSUBROUTINE