
      SUBROUTINE mass_diff_mod
!      
!     calculate mass diffusivity using Fuji et al.(1977) model 2015.07.27 JHLee (SNU)
!     calculate mass diffusivity using Chapman-Enskog theory
!
      USE VOL_DATA
      USE Zncg         , ONLY: ncg_diff,n_ncg_sp,ncg_species
      USE Zzone        , ONLY: ncell_fluid
!      
      IMPLICIT NONE
!      
      INTEGER i
      INTEGER ncg1,ncg2
!      
      REAL(8) mdiff,schmidt_turb,Sct
      REAL(8) Mols,sigmas,eks,Moln(8),sigman(8),ekn(8)
      REAL(8) c1,c2,c3,c4,c5,c6,c7,c8
      REAL(8) Molsr,Molnr1,Molnr2,Msn1,sigma1
      REAL(8) eksn1,Tstar1,omega1
      REAL(8) Dnum,Dden,D
!
!     NC gases are: Helium, Hydrogen, Nitrogen, Krypton, Xenon, Air, Argon, SF6. 
      DATA Mols/18.02d0/,sigmas/2.641d0/,eks/809.1d0/          ! for steam
      DATA Moln/4.002598d0, 2.01593d0, 28.01403d0, 83.800d0,   &
                131.300d0, 28.963d0, 39.948d0, 146.05d0/     ! molecular weight of NCG
      DATA sigman/2.551d0, 2.827d0, 3.798d0, 3.655d0,   &
                  4.047d0, 3.711d0, 3.542d0, 5.128d0/        ! characteristic length of NCG
      DATA ekn/10.22d0, 59.7d0, 71.4d0, 178.9d0,   &
               231.0d0, 78.6d0, 93.3d0, 222.1d0/             ! characteristic energy/Boltzmann's const. of NCG
      DATA c1,c2,c3,c4,c5,c6,c7,c8/1.06036d0,0.15610d0,0.19300d0,0.47635d0,   &
                                   1.03587d0,1.52996d0,1.76474d0,3.89411d0/
!
      IF(ncg_diff.eq.1)THEN
         DO i=1,ncell_fluid
            schmidt_turb=1.0d0 ! Schmidt number for turbulent diffusivity
            mdiff=7.65e-5*(cell%tg(i))**(11.0d0/6.0d0)/cell%p(i)+cell%tviscosg(i)/(cell%rhog(i)*schmidt_turb)
            cell%mdiff(i)=mdiff         
         ENDDO
!
      ELSEIF(ncg_diff.eq.2)THEN
         IF(n_ncg_sp.eq.1)THEN        ! binary D (steam-NCG)
            ncg1=ncg_species(1)
            DO i=1,ncell_fluid
               Molsr=1.0d0/Mols
               Molnr1=1.0d0/Moln(ncg1)
               Msn1=2.0d0/(Molsr+Molnr1)
!
               sigma1=(sigmas+sigman(ncg1))/2.0d0
!
               eksn1=DSQRT(eks*ekn(ncg1))
               Tstar1=cell%tg(i)/eksn1
!
               omega1=c1/Tstar1**c2
               omega1=omega1+c3/DEXP(c4*Tstar1)
               omega1=omega1+c5/DEXP(c6*Tstar1)
               omega1=omega1+c7/DEXP(c8*Tstar1)
!
               Dnum=0.00266d0*cell%tg(i)**1.5d0
               Dden=(cell%p(i)*1.0d-5)*DSQRT(Msn1)*omega1*sigma1*sigma1  ! unit conversion (Pa ---> bar)
               D=Dnum/Dden
               mdiff=D*1.0d-4     ! unit conversion (cm^2/s ---> m^2/s)
!
               Sct=0.7d0                ! turbulent schmidt number : 0.7(default value in Fluent 6.3)
               cell%mdiff(i)=mdiff+cell%tviscosg(i)/(cell%rhog(i)*Sct)
            ENDDO
         ELSEIF(n_ncg_sp.eq.2)THEN
            ncg1=ncg_species(1)
            ncg2=ncg_species(2)
            DO i=1,ncell_fluid
               IF(cell%quala(i).gt.0.99d0)THEN     ! binary D (NCG1-NCG2)
                  Molnr1=1.0d0/Moln(ncg1)
                  Molnr2=1.0d0/Moln(ncg2)
                  Msn1=2.0d0/(Molnr1+Molnr2)
!
                  sigma1=(sigman(ncg1)+sigman(ncg2))/2.0d0
!
                  eksn1=DSQRT(ekn(ncg1)*ekn(ncg2))
                  Tstar1=cell%tg(i)/eksn1
!
                  omega1=c1/Tstar1**c2
                  omega1=omega1+c3/DEXP(c4*Tstar1)
                  omega1=omega1+c5/DEXP(c6*Tstar1)
                  omega1=omega1+c7/DEXP(c8*Tstar1)
!
                  Dnum=0.00266d0*cell%tg(i)**1.5d0
                  Dden=(cell%p(i)*1.0d-5)*DSQRT(Msn1)*omega1*sigma1*sigma1  ! unit conversion (Pa ---> bar)
                  D=Dnum/Dden
                  mdiff=D*1.0d-4     ! unit conversion (cm^2/s ---> m^2/s)
               ENDIF
               Sct=0.7d0                ! turbulent schmidt number : 0.7(default value in Fluent 6.3)
               cell%mdiff(i)=mdiff+cell%tviscosg(i)/(cell%rhog(i)*Sct)
            ENDDO
         ENDIF
      ENDIF
!
      RETURN
      END SUBROUTINE mass_diff_mod
