!
      SUBROUTINE udfn_sg_outlet_property
!
      USE VOL_DATA                 
      USE Zbc_index    , ONLY: npb
      USE Zzone        , ONLY: ncell_fluid
!
      IMPLICIT NONE
!
      INTEGER i
!
!.....Set the void fraction at the steam separator to 1.0
!
      DO i=1,ncell_fluid
         IF(npb(i).eq.2)THEN
            cell%alphal(i)=0.0d0
            cell%alphag(i)=1.0d0
            cell%alphal_o(i)=0.0d0
            cell%alphag_o(i)=1.0d0
         ENDIF
      ENDDO
!
      RETURN
      END SUBROUTINE udfn_sg_outlet_property
