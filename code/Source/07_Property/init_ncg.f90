      SUBROUTINE init_ncg 
!
!     This routine initializes noncondensible gas DATA      
!                                                                       
      USE STM_TBL_cupid  , ONLY: wmole,dcvax,cvaox,uaox,rmolg, &
                                 tao,cvao,uao,dcva,ra
!
      IMPLICIT NONE
!
      INTEGER i      
!
!     REAL(8) tre(8),visa(8),thcax(8),  &
!             thcbx(8),advnon(8),dcon(8),qn(8)                                 
      REAL(8) qn(8)
!
!     wmolea  molecular mass of non-condensible gas.
!     rax     gas constant of non-condensible gas.
!     dcvax   same as above.
!     cvaox   same as above.
!     uaox    term in u = uao + integral (cv*dt) where u is internal
!             energy.
!     tao     term in cv = cvao + dcva*(t - tao) where cv is heat capacity
!             and t is temperature.
!     dconst  DIFfusion coefficient at reference conditions for non-
!             condensible gasses and steam.
!     noncn   number of non-condensible gasses.
!     prop    array for sth2x CALLs, also USEd for scratch.
!     s       same as above.
!
!     CHARACTER type(8)*8
!
!     DATA type/'helium','hydrogen','nitrogen','krypton','xenon','air', &
!               'argon','sf6'/
!     DATA tre/80.3d0,83.0d0,102.7d0,188.0d0,252.01d0,114.0d0,147.0d0,  &
!              0.0d0/                                                            
!     DATA visa/1.473d-6,6.675d-7,1.381d-6,2.386d-6,3.455d-6,1.492d-6,  &
!               1.935d-6,2.306654d-6/                                             
!     DATA thcax/2.639d-3,1.097d-3,5.314d-4,8.247d-5,4.351d-5,1.945d-4, &
!                2.986d-4,2.374568d-2/                                             
!     DATA thcbx/0.7085d0,0.8785d0,0.6898d0,0.8363d0,0.8616d0,0.8586d0, &
!                0.7224d0,0.0d0/                                                   
!     DATA advnon/2.67d0,6.12d0,18.5d0,24.5d0,32.7d0,19.7d0,16.2d0,     &
!                 71.3d0/
!     DATA dcon/3.9934d-4,4.2941d-4,1.2232d-4,9.5060d-5,8.2510d-5,      &
!               1.1886d-4,1.2039d-4,5.9826d-5/ 
!  
      DATA qn/0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,1.0d0,0.0d0,0.0d0/
!                                                                       
!     Data Statement:  Constant for evaluation of the DIFfusion            
!                      Coefficient of NC gas in Water Vapor.               
!                                                                         
!                      diffc = Dconst * T**1.75 / P                        
!                                                                       
!     Ref:  eq. 11-4.1 of "Properties of Gases and Liquids"           
!           by Reid, Praudnitz & Sherwood.                            
!           3rd edition, McGraw-Hill Book Co., 1977.                  
!                                                                       
!     NC gases are: Helium, Hydrogen, Nitrogen, Krypton, Xenon,         
!                   Air, Argon, SF6.                                    
!                                                                       
      tao=250.0d0 
      cvao=0.0d0 
      uao=0.0d0 
      dcva=0.0d0 
      ra=0.0d0 
!
      DO i=1,8 
         cvao=cvao+cvaox(i)*qn(i)
         uao=uao+uaox(i)*qn(i)
         dcva=dcva+dcvax(i)*qn(i)
         ra=ra+rmolg/wmole(i)*qn(i)
      END DO 
!
      END SUBROUTINE init_ncg 




