!
      SUBROUTINE initialize_least_square_option
!
!     Define a cell on walls to avoid least square method
!    
      USE Zmpi         , ONLY: maxmt_fluid
      USE Zzone        , ONLY: ncell_fluid
      USE Zparam       , ONLY: ndim
      USE Znum_cell    , ONLY: i_neigh
      USE Zbc_index    , ONLY: nbcon
      USE Zgrad_ls_c3d , ONLY: lsindex
      USE Zvec_geo     , ONLY: xn_nf
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER i,j,wallopt
      INTEGER imsi
!.....Local arrays
      REAL(8):: xn(maxmt_fluid)
!            
!
!.....Find bottom axis with gravity vector
! 
!     DO ix=1,ndim
!        IF(grav(ix)*grav(ix).gt.0.5d0)EXIT
!     ENDDO      
!      
      IF(ndim.eq.2) THEN
         CALL get_scalar_variable_n(xn_nf(1,2),xn) 
         DO i=1,ncell_fluid
            wallopt=0
            imsi=1
            DO j=i_neigh(i),i_neigh(i+1)-1
               IF(nbcon(j).eq.-1)THEN
                  wallopt=1
                  IF(abs(xn(j)).gt.1.0d-3)THEN 
                     imsi=0                       !bottom or top wall
                     GOTO 100
                  ENDIF
               ENDIF             
            ENDDO
100         CONTINUE
            IF(wallopt.eq.1)THEN
               IF(imsi.eq.0)THEN
!
!................Use Least Square            
!
                  lsindex(i)=1
              ELSEIF(imsi.eq.1)THEN
!
!................Avoid Least Square            
!         
                  lsindex(i)=0
              ENDIF   
            ENDIF   
         ENDDO  
      ELSE
         CALL get_scalar_variable_n(xn_nf(1,3),xn) 
         DO i=1,ncell_fluid
            wallopt=0
            imsi=1
            DO j=i_neigh(i),i_neigh(i+1)-1
               IF(nbcon(j).eq.-1)THEN
                  wallopt=1
                  IF(abs(xn(j)).gt.1.0d-3)THEN 
                     imsi=0                       !bottom or top wall
                     GOTO 200
                  ENDIF
               ENDIF             
            ENDDO
200         CONTINUE
!         
            IF(wallopt.eq.1)THEN
               IF(imsi.eq.0)THEN
!
!...........Use Least Square            
!
                  lsindex(i)=1
              ELSEIF(imsi.eq.1)THEN
!
!...........Avoid Least Square            
!         
                  lsindex(i)=0
              ENDIF   
            ENDIF   
         ENDDO  
      ENDIF
!
      END SUBROUTINE initialize_least_square_option
