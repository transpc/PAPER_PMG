!
      SUBROUTINE rv_int_ht
!
      USE VOL_DATA       , ONLY: cell
      USE Zzone          , ONLY: ncell_fluid         
!
      IMPLICIT NONE
!
      INTEGER :: i,ag_option
      REAL(8),DIMENSION(ncell_fluid) :: agi
!
!.....initialization
!
      ag_option=0       
!
!.....void fraction definition
!
!!!!!removed no need
!     IF(ag_option.eq.1) CALL alpha_predictor(ag,al)
      DO i=1,ncell_fluid
         agi(i)=cell%alphag(i)
      ENDDO
!
      CALL rv_ihtc_main(agi)
!      
      END SUBROUTINE rv_int_ht
