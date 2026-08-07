!
      SUBROUTINE wall_condensation_multi_diff
!
!     This routine calculates binary diffusivity.
!
      USE VOL_DATA
      USE Zbc_index       , ONLY: icell_type
      USE Zmodel          , ONLY: coef_diff
      USE Zncg            , ONLY: advn_cell,wmole_gas
      USE Zzone           , ONLY: ncell_fluid
!
      IMPLICIT NONE
!
      INTEGER i
!
      REAL(8) molh2o,advh2o,dcnst
  
!
      DATA molh2o,advh2o/18.0d0,12.7d0/       ! molecular weights of water vapor & air
!
      coef_diff(:)=0.0d0
!                                                                       
!.....Mass diffusivity is calculated using eqn.11-4.4 (Fuller) in the properties of gases and liquid by reid,praudnitz,sherwood,third ed. mc-graw-hill book co.,1977.                                
!       
      DO i=1,ncell_fluid
         IF(icell_type(i).eq.1)THEN      
            dcnst=(sqrt(1.0d0/molh2o+1.0d0/wmole_gas(i)))/((advh2o**0.33333d0+advn_cell(i)**0.33333d0)**2)                                                    
            coef_diff(i)=0.0101325d0*dcnst*cell%tg(i)**1.75d0/cell%p(i)     
         ENDIF
      ENDDO   
!
      RETURN
      END SUBROUTINE wall_condensation_multi_diff   

!!
!      SUBROUTINE wall_condensation_binary_diff
!!
!!     This routine calculates binary diffusivity.
!!
!      USE VOL_DATA
!      USE Zbc_index       , ONLY: icell_type
!      USE Zmodel          , ONLY: binaryD
!      USE Zzone           , ONLY: ncell_fluid
!!
!      IMPLICIT NONE
!!
!      INTEGER i
!!
!      LOGICAL D_STP,D_CE
!!
!      REAL(8) T0,P0,TgT0,P0P
!      REAL(8) Molv,Mola,Molvr,Molar,Mva
!      REAL(8) sigma,sigmav,sigmaa
!      REAL(8) ekva,ekv,eka,Tstar
!      REAL(8) omega,c1,c2,c3,c4,c5,c6,c7,c8,D,D1,D2      
!!
!      DATA Molv,Mola/18.02d0,28.97d0/       ! molecular weights of water vapor & air
!      DATA sigmav,sigmaa/2.641d0,3.711d0/   ! characteristic length of water vapor & air
!      DATA ekv,eka/809.1d0,78.6d0/          ! characteristic energy/Boltzmann's const. of water vapor & air
!      DATA c1,c2,c3,c4,c5,c6,c7,c8/1.06036d0,0.15610d0,0.19300d0,0.47635d0,   &
!                                   1.03587d0,1.52996d0,1.76474d0,3.89411d0/
!!
!      D_STP=.true.
!      D_CE=.false.
!!
!      binaryD(:)=0.0d0
!!
!      DO i=1,ncell_fluid
!         IF(icell_type(i).eq.1)THEN
!            IF(D_STP)THEN         ! bindary diffusivity at near Standard Temp. and Pressure
!               T0=298.0d0
!               P0=1.0d5
!               TgT0=cell%tg(i)/T0
!               P0P=P0/cell%p(i)
!               binaryD(i)=(2.6d0*1.0d-5)*(TgT0**1.75)*(P0P)
!            ELSEIF(D_CE)THEN       ! Chapman-Enskog theory
!               Molvr=1.0d0/Molv
!               Molar=1.0d0/Mola
!               Mva=2.0d0/(Molvr+Molar)
!!
!               sigma=(sigmav+sigmaa)/2.0d0
!!
!               ekva=DSQRT(ekv*eka)
!               Tstar=cell%tg(i)/ekva
!!
!               omega=c1/Tstar**c2
!               omega=omega+c3/DEXP(c4*Tstar)
!               omega=omega+c5/DEXP(c6*Tstar)
!               omega=omega+c7/DEXP(c8*Tstar)
!!
!               D1=0.00266d0*cell%tg(i)**1.5d0
!               D2=(cell%p(i)*1.0d-5)*DSQRT(Mva)*omega*sigma*sigma  ! unit conversion (Pa ---> bar)
!               D=D1/D2
!               binaryD(i)=D*1.0d-4     ! unit conversion (cm^2/s ---> m^2/s)
!            ENDIF
!         ENDIF
!      ENDDO
!!
!      RETURN
!      END SUBROUTINE wall_condensation_binary_diff
