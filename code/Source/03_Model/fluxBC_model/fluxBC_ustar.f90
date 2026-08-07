!
      SUBROUTINE fluxBC_ustar
!
! update u* when choking happens      
!
      USE VOL_DATA     
      USE Zparam       , ONLY: ndim
      USE Zvector      , ONLY: vl_n,vg_n
      USE Zpress       , ONLY: p,dpdx
      USE Zpress_coeff , ONLY: coefp_g,coefp_l
      USE Zgradoption  , ONLY: grav_grad
      USE Zcore        , ONLY: np
      USE Zmpi         , ONLY: ncell_fp
      USE Zzone        , ONLY: ncell_fluid
      USE Zconst2      , ONLY: dt
      USE Zrv_choke    , ONLY: time_cflow_on,choke_update,choke
      USE Ztimecon     , ONLY: time
      USE Zvalve       , ONLY: time_valve_closed
      USE Zrv_model    , ONLY: rv_mcp,rv_valve
!
      IMPLICIT NONE
!      
      INTEGER :: i
      REAL(8) :: time_valve_max
      REAL(8) :: ptmp(ncell_fp),dppdx(ncell_fp,ndim)
!      
!......Apply Neumann BC         
!      
!      check flow is choked
      choke_update=0
      IF(choke) choke_update=1
      IF(np.gt.1) CALL allreducei_i1(choke_update)       
!
!      check valve is closed.
      time_valve_max=MAXVAL(time_valve_closed(:,2))
      IF(rv_valve.eq.0) time_valve_max=0.d0
!      
!      check mcp is equiped: rv_mcp=1
!      (blank)
!      
!......Apply Neumann BC and Update U*        
!       
      IF((choke_update.eq.1.and.time.ge.time_cflow_on).or. &      !choked condition
         (time.le.time_valve_max).or. & !valve model is applied 
         (rv_mcp.eq.1) ) THEN          !mcp model is applied
!      
!.........Apply Neumann BC for the flux-controlled faces.
!           
         ptmp=p
         IF(grav_grad.eq.0)THEN
            CALL grad_press(ptmp,dppdx,0)
         ELSEIF(grav_grad.eq.1)THEN
            CALL grad_pressK1(ptmp,dppdx,0)
         ELSEIF(grav_grad.ge.2)THEN
             Print*,'>>> Error'
             Print*,'>>> grav_grad for choke model is 0 or 1'
             PAUSE
             STOP
         ENDIF      
!      
!.........Update U*          
!   
         IF(ndim.eq.2) THEN      
            DO i=1,ncell_fluid
               vl_n(i,1)=vl_n(i,1)+coefp_l(i)*(dpdx(i,1)-dppdx(i,1))*dt
               vl_n(i,2)=vl_n(i,2)+coefp_l(i)*(dpdx(i,2)-dppdx(i,2))*dt
               vg_n(i,1)=vg_n(i,1)+coefp_g(i)*(dpdx(i,1)-dppdx(i,1))*dt
               vg_n(i,2)=vg_n(i,2)+coefp_g(i)*(dpdx(i,2)-dppdx(i,2))*dt
            ENDDO
         ELSE
            DO i=1,ncell_fluid
               vl_n(i,1)=vl_n(i,1)+coefp_l(i)*(dpdx(i,1)-dppdx(i,1))*dt
               vl_n(i,2)=vl_n(i,2)+coefp_l(i)*(dpdx(i,2)-dppdx(i,2))*dt
               vl_n(i,3)=vl_n(i,3)+coefp_l(i)*(dpdx(i,3)-dppdx(i,3))*dt
               vg_n(i,1)=vg_n(i,1)+coefp_g(i)*(dpdx(i,1)-dppdx(i,1))*dt
               vg_n(i,2)=vg_n(i,2)+coefp_g(i)*(dpdx(i,2)-dppdx(i,2))*dt
               vg_n(i,3)=vg_n(i,3)+coefp_g(i)*(dpdx(i,3)-dppdx(i,3))*dt             
            ENDDO          
         ENDIF
!         
      ENDIF
!     
      RETURN
      END SUBROUTINE fluxBC_ustar

    