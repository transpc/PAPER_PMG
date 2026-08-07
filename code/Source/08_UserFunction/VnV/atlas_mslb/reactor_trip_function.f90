      Function reactor_trip_function(rx_trip_signal)
!      
!      USE Zmars    ,ONLY: rx_trip_mars
!      
      IMPLICIT NONE
      INCLUDE '../../../10_LinkToMARS/c3com.h' !rx_trip_mars
      
!
      INTEGER reactor_trip_function,rx_trip_signal
!      
      IF(rx_trip_signal.eq.1)THEN
         reactor_trip_function=1
         RETURN      
      ENDIF   
!
      IF(i3rx_trip.gt.0)then
         reactor_trip_function=1
      ELSE
         reactor_trip_function=0
      ENDIF
!       
      END FUNCTION reactor_trip_function      
