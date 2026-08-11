!
!     communicate_allb — 초기화 시 주요 변수 일괄 고스트 교환 (C010-1 에서 communicate.f90 로부터
!     verbatim 분리: CUPID 본체 모듈(VOL_DATA/Zpress/Zvector/Zare) 의존이 이 루틴뿐이라
!     PMG 단독 하네스가 communicate.f90 원본을 그대로 컴파일할 수 있도록 파일 단위를 나눔)
!
      SUBROUTINE communicate_allb
!DEC$IF defined (mpi_flag)
!
!     Comunicate major variables for parallel computing
!
      USE Zinterface
      USE VOL_DATA , ONLY: cell
      USE Zmpi     , ONLY: ncell_fp
      USE Zpress   , ONLY: p
      USE Zvector  , ONLY: vg_o,vl_o,vd_o
      USE Zare     , ONLY: are_gas,are_liq
!
      IMPLICIT NONE
!
!.....Local variables
      LOGICAL,SAVE :: initial=.TRUE.
      INTEGER :: i
!
      IF(initial)THEN
         CALL communicate_1d(are_liq, &
                             are_gas, &
                             cell%eg)
         initial=.false.
      ENDIF
      CALL communicate_1d(cell%alphag, &
                          cell%alphal, &
                          cell%alphad, &
                          cell%quala,  &
                          cell%tl_o,   &
                          cell%tg_o,   &
                          cell%condl,  &
                          cell%condg)
      CALL communicate_1d(cell%rhom,  &
                          cell%rhomr, &
                          cell%el,    &
                          cell%hg,    &
                          cell%hl,    &
                          cell%rhog,  &
                          cell%rhol,  &
                          cell%ha)
      CALL communicate_1d(p)
      CALL communicate_2d(vl_o, &
                          vg_o, &
                          vd_o)
!
!.....For mass_conv, erg_conv
!
      DO i=1,ncell_fp
         cell%alphag_o(i)=cell%alphag(i)
         cell%alphal_o(i)=cell%alphal(i)
         cell%alphad_o(i)=cell%alphad(i)
      ENDDO
!DEC$ELSE
!
      USE VOL_DATA , ONLY: cell
      USE Zzone    , ONLY: ncell_fluid
!
      IMPLICIT NONE
!.....Local variables
      INTEGER :: i
!
!.....For mass_conv, erg_conv
!
      DO i=1,ncell_fluid
         cell%alphag_o(i)=cell%alphag(i)
         cell%alphal_o(i)=cell%alphal(i)
         cell%alphad_o(i)=cell%alphad(i)
      ENDDO
!DEC$ENDIF
!
      END SUBROUTINE communicate_allb
