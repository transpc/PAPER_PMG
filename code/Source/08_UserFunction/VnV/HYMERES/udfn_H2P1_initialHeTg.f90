!
      SUBROUTINE udfn_H2P1_initialHeTg
!
      USE VOL_DATA
      USE Zconst1      , ONLY: vv_prob
      USE Zconst2      , ONLY: grav
      USE Zcoord1      , ONLY: xloc
      USE Zzone        , ONLY: ncell_fluid
!
      IMPLICIT NONE
!
      INTEGER :: i,j
      REAL(8) :: yi_he,yi_steam,mol_he,mol_steam,mi_he,mi_steam,xi_he,FACT
      REAL(8),ALLOCATABLE::initialHe(:,:),initialTg(:,:)
!
      ALLOCATE(initialHe(2,19),initialTg(2,16))
      IF(vv_prob.eq.'VD_h2p1_0')THEN
         OPEN(873,file='initialHe.dat',status='old',form='unformatted')         !initial He
         DO i=1,19
            READ(873) initialHe(1,i), initialHe(2,i)
         ENDDO
      ELSE
         OPEN(873,file='initialHe.dat',status='old')         !initial He
         READ(873,*) initialHe(:,:)
      ENDIF
      CLOSE(873)
!
      mol_he=4.003d0
      mol_steam=18.02d0
      DO i=1,ncell_fluid
         IF(grav(2).gt.grav(3))THEN         !height:z
            IF(xloc(i,3).le.initialHe(1,1))THEN
               yi_he=initialHe(2,1)
            ELSE
               DO j=1,18
                  IF(xloc(i,3).gt.initialHe(1,j).and.xloc(i,3).le.initialHe(1,j+1))THEN
                     FACT=(xloc(i,3)-initialHe(1,j))/(initialHe(1,j+1)-initialHe(1,j))
                     yi_he=initialHe(2,j)+FACT*(initialHe(2,j+1)-initialHe(2,j))
                  ENDIF
               ENDDO
            ENDIF
         ELSEIF(grav(2).lt.grav(3))THEN     !height:y
            IF(xloc(i,2).le.initialHe(1,1))THEN
               yi_he=initialHe(2,1)
            ELSE
               DO j=1,18
                  IF(xloc(i,2).gt.initialHe(1,j).and.xloc(i,2).le.initialHe(1,j+1))THEN
                     FACT=(xloc(i,2)-initialHe(1,j))/(initialHe(1,j+1)-initialHe(1,j))
                     yi_he=initialHe(2,j)+FACT*(initialHe(2,j+1)-initialHe(2,j))
                  ENDIF
               ENDDO
            ENDIF
         ENDIF
         yi_steam=1.0d0-yi_he
         mi_he=yi_he*mol_he
         mi_steam=yi_steam*mol_steam
         xi_he=mi_he/(mi_he+mi_steam)
         cell%quala(i)=xi_he
      ENDDO
!
      IF(vv_prob.eq.'VD_h2p1_0')THEN
         OPEN(874,file='initialTg.dat',status='old',form='unformatted')         !initial Tg
         DO i=1,16
            READ(874) initialTg(1,i), initialTg(2,i)
         ENDDO
      ELSE
         OPEN(874,file='initialTg.dat',status='old')         !initial Tg
         READ(874,*) initialTg(:,:)
      ENDIF
      CLOSE(874)
!
      DO i=1,ncell_fluid
         IF(grav(2).gt.grav(3))THEN         !height:z
            IF(xloc(i,3).le.initialTg(1,1))THEN
               cell%tg(i)=initialTg(2,1)
            ELSE
               DO j=1,15
                  IF(xloc(i,3).gt.initialTg(1,j).and.xloc(i,3).le.initialTg(1,j+1))THEN
                     FACT=(xloc(i,3)-initialTg(1,j))/(initialTg(1,j+1)-initialTg(1,j))
                     cell%tg(i)=initialTg(2,j)+FACT*(initialTg(2,j+1)-initialTg(2,j))
                  ENDIF
               ENDDO
            ENDIF
         ELSEIF(grav(2).lt.grav(3))THEN     !height:y
            IF(xloc(i,2).le.initialTg(1,1))THEN
               cell%tg(i)=initialTg(2,1)
            ELSE
               DO j=1,15
                  IF(xloc(i,2).gt.initialTg(1,j).and.xloc(i,2).le.initialTg(1,j+1))THEN
                     FACT=(xloc(i,2)-initialTg(1,j))/(initialTg(1,j+1)-initialTg(1,j))
                     cell%tg(i)=initialTg(2,j)+FACT*(initialTg(2,j+1)-initialTg(2,j))
                  ENDIF
               ENDDO
            ENDIF
         ENDIF
         cell%tg(i)=cell%tg(i)+273.16d0
      ENDDO
!
      DEALLOCATE(initialHe,initialTg)
!
      RETURN
      END SUBROUTINE udfn_H2P1_initialHeTg
!
!------------------------------------------------------------------------
! hymeres2
!
! <input>
! 870 : injectionT.dat
! 871 : injectionV.dat
! 872 : injectionTKE.dat
! 873 : initialHe.dat
! 874 : initialTg.dat
! 875 : inputTw.dat
! 876 : inputTp.dat
!
! <output>
! 442 : z_He_initial.dat
! 443 : z_Tg_initial.dat
! 444 : z_He_transient20.dat
! 445 : z_He_transient14.dat
! 446 : z_He_transient26.dat
! 447 : z_Tg_transient20.dat
! 448 : z_Tg_transient14.dat
! 449 : z_Tg_transient26.dat
! 451 : z_Tg_erosionfront.dat
! 471-479 : FOV
! 571-579 : TE
! 671-679 : Vy
! 771-779 : avg time FOV
! 871-879 : avg time Vy

