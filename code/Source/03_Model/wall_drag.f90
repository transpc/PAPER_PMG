!
      SUBROUTINE wall_drag
!     
!     This calculates wall drag coefficient in porous cells
!
      USE VOL_DATA     , ONLY: cell
      USE Zmpi         , ONLY: ncell_fp
      USE Zzone        , ONLY: nzone,ncell_fluid
      USE Zconst2      , ONLY: hydraulicd
      USE Zporous      , ONLY: udfi_subchannel_flowdir,fric_model_gas,fric_model_liq
      USE Zmodel       , ONLY: s_wall_fric
      USE Zporous      , ONLY: s_subchannel_fric
      USE Zvector      , ONLY: ug_o,ul_o
!      USE Zcoord3      , ONLY: porosity
!
      IMPLICIT NONE
!      
      INTEGER i
!      
      REAL(8) vfwg_i,vfwl_i
      REAL(8) fric,Rel,Reg
      REAL(8) Dh,Re,fwk1,fwk2,fwk3
!     local arrays
      REAL(8) :: vfwg(ncell_fp),vfwl(ncell_fp)
!
      IF(s_wall_fric.eq.'none')RETURN       
!
      IF(s_wall_fric.eq.'udf')THEN            
         CALL udfn_wall_drag(vfwl,vfwg) 
      ENDIF
      DO i=1,ncell_fluid
!      
!........Default value is zero      
!
         vfwg_i=0.0d0
         vfwl_i=0.0d0
         dh=hydraulicd(i)
!          
!........Use of user defined function
!
         IF(s_wall_fric.eq.'udf')THEN            
!           CALL udfn_wall_drag(i,vfwl_i,vfwg_i) 
            vfwg_i=vfwg(i)
            vfwl_i=vfwl(i)
!      
!........Kakac model      
!
         ELSEIF(s_wall_fric.eq.'kakac') THEN
            IF(cell%regime(i).eq.11 .or. cell%regime(i).eq.12 .or. cell%regime(i).eq.21 .or. &
            cell%regime(i).eq.22 .or. cell%regime(i).eq.3) THEN
               Rel=DMAX1(1.0d0,(cell%rhol(i)*ul_o(i)*(2.0d0*HydraulicD(i))/cell%lviscosl(i)))
               CALL mom_wall_kakac(Rel,fric)  
               vfwl_i=fric/(2.0d0*0.11d0)*cell%rhol(i)*ul_o(i)
            ELSEIF(cell%regime(i).eq.13)THEN
               Reg=DMAX1(1.0d0,cell%rhog(i)*ug_o(i)*(2.0d0*HydraulicD(i))/cell%lviscosg(i))
               CALL mom_wall_kakac(Reg,fric)             
               vfwg_i=fric/(2.0d0*0.11d0)*cell%rhog(i)*ug_o(i)
               vfwl_i=cell%alphal(i)*fric/(2.0d0*HydraulicD(i))*cell%rhol(i)*ul_o(i)
            ENDIF
!      
!........Darcy's formation   
!   
         ELSEIF(s_wall_fric.eq.'darcy') THEN 
!     
!...........For liquid phase  
!
            Re=DMAX1(50.d0,cell%alphal(i)*cell%rhol(i)*ul_o(i)*Dh/cell%lviscosl(i))
            fwk1=16.d0/Re*4.d0
            fwk2=0.0791d0/(Re**0.25d0)*4.d0
            fwk3=0.0008d0+0.05525d0/(Re**0.237d0)*4.d0
!
!...........Sorting maximum fwki(i=1,2,3)      
!            
            IF(fwk1.gt.fwk2) THEN        
               cell%fwkl(i)=fwk1
            ELSE
               cell%fwkl(i)=fwk2
            ENDIF
            IF(cell%fwkl(i).le.fwk3) THEN
               cell%fwkl(i)=fwk3
            ENDIF
!            
            vfwl_i=cell%fwkl(i)/2.d0/Dh*cell%rhol(i)*ul_o(i) * cell%alphal(i)**(1.d0/8.d0) 
!     
!...........For gas phase
!
            Re=DMAX1(50.d0,cell%alphag(i)*cell%rhog(i)*ug_o(i)*Dh/cell%lviscosg(i))
            fwk1=16.d0/Re*4.d0
            fwk2=0.0791d0/(Re**0.25d0)*4.d0
            fwk3=0.0008d0+0.05525d0/(Re**0.237d0)*4.d0
!
!...........Sorting maximum fwki(i=1,2,3)      
!        
            IF(fwk1.gt.fwk2) THEN
               cell%fwkg(i)=fwk1
            ELSE
               cell%fwkg(i)=fwk2
            ENDIF
            IF(cell%fwkg(i).le.fwk3) THEN
               cell%fwkg(i)=fwk3
            ENDIF
!            
            vfwg_i=cell%fwkg(i)/2.d0/Dh*cell%rhog(i)*ug_o(i) * cell%alphag(i)**(1.d0/8.d0)        
            IF(nzone(i).ge.2)  THEN
               vfwl_i= vfwl_i*50.d0
               vfwg_i= vfwg_i*50.d0
            ENDIF 
!      
!........isotropic friction loss part of Anisotropic Friction Loss Model
!
         ELSEIF(s_subchannel_fric.eq.'aniso_fric_imp') THEN
            vfwl_i=-fric_model_liq(i,udfi_subchannel_flowdir)
            vfwg_i=-fric_model_gas(i,udfi_subchannel_flowdir)
         ELSE
            vfwg_i=0.0d0
            vfwl_i=0.0d0          
         ENDIF
!
      ENDDO
!     
      END SUBROUTINE wall_drag
