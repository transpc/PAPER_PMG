      SUBROUTINE rv_ruplas(sigma,dtdt,tplask,trupk) 
!
!deck ruplas                                                            
!                                                                       
!   calculate plastic deformation and rupture temperatures              
!
      IMPLICIT none 
!                                                                       
      REAL(8) sigma,dtdt,tplask,trupk 
!                                                                       
! ********* input *************                                         
!     sigma  = hoop stress (pa)                                        
!     dtdt   = heatup rate (k/s)                                       
! ******** output ***************                                       
!     trupk  = rupture temperature (k)                                 
!     tplask = plastic deformation temperature (k)                     
!                                                                       
!     rupture temperature is from equation 3-2, page 10 in nureg-0630     
!                                                                       
!.....Local variables.                                                     
!
      REAL(8) degkc,delt,heatup,pakpsi,rate,sigmak,tplas,trup 
!
!.....Conversion from pa to kpsi                 
!
      DATA pakpsi/1.45037747d-7/ 
!      
!.....Conversion from K to oC                     
!
      DATA degkc/273.2d0/ 
!                            
!.....Base heat-up rate deg oC/s
!
      DATA heatup/28.0d0/ 
!      
!.....Convert sigma from pa to kpsi
!
      sigmak=sigma*pakpsi 
!
!.....Calc. heatup rate, ratio to 28 k/sec (limit to range between 0 and 1)
!
      rate=DMIN1(1.0d0,DMAX1(0.0d0,dtdt/heatup)) 
!      
!.....Compute rupture temperature (oC)            
!
      trup=3960.d0-20.4d0*sigmak/(1.0d0+rate)-8510000.0d0*sigmak/(100.0d0*(1.0d0+rate)+2790.0d0*sigmak)   !Tf: Rupture Temperature
!                                                                       
!.....Check that ballooning does not occur during compressive stress                  
!
      IF(sigma.lt.0.0d0) trup=1300.0d0 
!                                                                       
!.....Calculate delta temperature (oC)            
!
      IF(trup.lt.700.0d0) THEN
         delt=70.0d0
      ELSEIF(trup.gt.1300.0d0) THEN
         delt=155.0d0
      ELSE
         delt=70.0d0+(trup-700.0d0)*(85.0d0/600.0d0)
      ENDIF
!                                                                       
!.....Find cladding temperature at which plastic deformation begins (oC)            
!
      tplas=trup-delt     !Tplas: Rupture Temperature
!      
!.....Convert to degrees (K)
!
      tplask=tplas+degkc 
      trupk=trup+degkc 
!      
      RETURN 
      END SUBROUTINE rv_ruplas   
!