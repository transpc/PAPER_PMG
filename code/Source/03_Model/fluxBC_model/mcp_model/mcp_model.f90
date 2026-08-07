!
      SUBROUTINE mcp_model(vol_flow)
!
!     MCP model with homologous curve
!
      USE Vol_DATA      ,ONLY: cell
      USE Zcore         ,ONLY: myrank,np
      USE Zpress        ,ONLY: p
      USE Zconst2       ,ONLY: dt
      USE Ztimecon      ,ONLY: time
      USE Zmcp          
      
      IMPLICIT NONE
     
      INTEGER :: zz,i,ii,kk,tt,j
      INTEGER,SAVE:: icnt_time
      REAL(8) :: agr,alr,cf,dz_duct,flux_err,hratio,spdrat,wratio,delt
      REAL(8),SAVE :: write_interval
      REAL(8) :: vol_flow(num_mcp)
      REAL(8),DIMENSION(:),SAVE,ALLOCATABLE:: agb,alb,rgb,rlb,p_mcp_in,p_mcp_out,dp_loop,rho_avg, &   !local
                                              agb_global,alb_global,rgb_global,rlb_global,p_mcp_in_global,p_mcp_out_global,dp_loop_global,rho_avg_global,rpm_controlled,head_pump,dp_err(:) !global
      LOGICAL,SAVE :: initial,initial_write
      DATA initial,initial_write /.true.,.true./
!
!......Sum of properties at MCP faces
!      
!      Initialize the properties       
      IF(initial) THEN
          zz=num_mcp
          ALLOCATE(agb(zz),alb(zz),rgb(zz),rlb(zz))
          ALLOCATE(p_mcp_in(zz),p_mcp_out(zz)) !1=MCP cells, 2=Distributor cells
          ALLOCATE(dp_loop(zz),rho_avg(zz)) !1=MCP cells, 2=Distributor cells
!         
          zz=num_mcp
          ALLOCATE(agb_global(zz),alb_global(zz),rgb_global(zz),rlb_global(zz))
          ALLOCATE(p_mcp_in_global(zz),p_mcp_out_global(zz)) !1=MCP cells, 2=Distributor cells
          ALLOCATE(dp_loop_global(zz),rho_avg_global(zz)) !1=MCP cells, 2=Distributor cells
!          
          ALLOCATE(head_mult(zz),dp_pump(zz),rated_vol_flow(zz),rpm_controlled(zz),head_pump(zz),dp_err(zz))
          initial=.false.
          
          head_mult=0.d0
          head_pump=0.d0
          dp_err=0.d0
      ENDIF
      agb=0.d0
      alb=0.d0
      rgb=0.d0
      rlb=0.d0
      p_mcp_in=0.d0
      p_mcp_out=0.d0
      dp_loop=0.d0
      rho_avg=0.d0
!
      agb_global=0.d0
      alb_global=0.d0
      rgb_global=0.d0
      rlb_global=0.d0
      p_mcp_in_global=0.d0
      p_mcp_out_global=0.d0
      dp_loop_global=0.d0
      rho_avg_global=0.d0      
!      
!      Sum of the properties      
      DO zz=1,num_mcploc
         tt=mapping_mcp(zz)
         DO i=1,num_mcpface(zz)
            ii=icell_mcp(zz,i)                 !MCP cell
            kk=ocell_mcp(zz,i)                 !Distributor cell 
            agr=cell%alphag(ii)*cell%rhog(ii)  
            alr=cell%alphal(ii)*cell%rhol(ii)  
            agb(tt)=agb(tt)+cell%alphag(ii)    !ag
            alb(tt)=alb(tt)+cell%alphal(ii)    !al  
            rgb(tt)=rgb(tt)+agr                !ag*rhog  
            rlb(tt)=rlb(tt)+alr                !al*rhol   
            p_mcp_in(tt)=p_mcp_in(tt)+p(ii)    !p (1=MCP,num_mcploc)
            p_mcp_out(tt)=p_mcp_out(tt)+p(kk)  !p (2=Distributor,num_mcploc) p(1:ncell_fp) <<- kk is available
         ENDDO
      ENDDO
      IF(np.gt.1) THEN
          CALL allreduce_r(agb,agb_global,num_mcp)
          CALL allreduce_r(alb,alb_global,num_mcp)
          CALL allreduce_r(rgb,rgb_global,num_mcp)
          CALL allreduce_r(rlb,rlb_global,num_mcp)
          DO i=1,num_mcp
             CALL allreduce_r(p_mcp_in(i),p_mcp_in_global(i),num_mcp)
             CALL allreduce_r(p_mcp_out(i),p_mcp_out_global(i),num_mcp)
          ENDDO   
      ELSE
          agb_global=agb
          alb_global=alb
          rgb_global=rgb
          rlb_global=rlb
          p_mcp_in_global=p_mcp_in
          p_mcp_out_global=p_mcp_out
      ENDIF
!
!      Averaging the properties (global domain)
      DO tt=1,num_mcp
         IF(alb_global(tt).gt.1.0d-8)THEN
            rlb_global(tt)=rlb_global(tt)/alb_global(tt) !void fraction average (alb)
         ENDIF
         IF(agb_global(tt).gt.1.0d-8)THEN
            rgb_global(tt)=rgb_global(tt)/agb_global(tt)
         ENDIF
         agb_global(tt)=agb_global(tt)/num_mcpface_global(tt)
         agb_global(tt)=DMIN1(DMAX1(0.0d0,agb_global(tt)),1.0d0)
         alb_global(tt)=1.0d0-agb_global(tt)
         p_mcp_in_global(tt)=p_mcp_in_global(tt)/num_mcpface_global(tt)
         p_mcp_out_global(tt)=p_mcp_out_global(tt)/num_mcpface_global(tt)
!         
         dp_loop_global(tt)=p_mcp_in_global(tt)-p_mcp_out_global(tt)
         rho_avg_global(tt)=agb_global(tt)*rgb_global(tt)+alb_global(tt)*rlb_global(tt)
      ENDDO     
!
!......MCP model      
!     
!      define steady or transient condition
      init_mcp=0 !transient
      IF(time.lt.mcp_transient_start(1)) init_mcp=1 !steady state
!
!      pump speed control
      rpm_controlled(:)=0.d0  !array=num_mcp (global domain)
      IF(init_mcp.eq.0) THEN
         DO i=1,num_mcp_transient
            IF(i.eq.num_mcp_transient) EXIT
            IF(time.ge.mcp_transient_start(i).and.time.lt.mcp_transient_start(i+1)) THEN
                delt=mcp_transient_start(i+1)-mcp_transient_start(i)
                DO j=1,num_mcp
                   rpm_controlled(j)=speed_pump(j,i)+(speed_pump(j,i+1)-speed_pump(j,i))/delt*(time-mcp_transient_start(i))
                ENDDO   
            ENDIF 
         ENDDO   
      ENDIF
!
!      MCP model      
      DO tt=1,num_mcp
         dz_duct=mcp_vol_global(tt)/mcp_area_global(tt)
         cf=dt*mcp_area_global(tt)/rho_avg_global(tt)/dz_duct    ! (s*m2*m3)/(kg*m)   !based on one face area
!
         IF(init_mcp.eq.1)THEN
            head_mult(tt)=-dp_loop_global(tt)/rho_avg_global(tt)/(9.8d0*rated_pump_hd(tt))  !averaged
            dp_pump(tt)=-dp_loop_global(tt)
            vol_flow(tt)=rwinit(tt)/rho_avg_global(tt)           !volumetric flow rate (m3/s) per one MCP
            rated_vol_flow(tt)=vol_flow(tt)                      !volumetric flow rate (m3/s) per one MCP  
         ELSE
!            mcp_stop=ALL(speed_pump(1:num_mcp).lt.1.0d0)
            spdrat=rpm_controlled(tt)/rated_pump_speed(tt)
            wratio=vol_flow(tt)/rated_vol_flow(tt)               !vol_flow and rated_vol_flow should be based on the same concept (one face or total face)
            CALL mcp_head(hratio,spdrat,wratio)                  !mcp head ratio 
!               
            head_pump(tt)=head_mult(tt)*hratio*rated_pump_hd(tt)
            dp_pump(tt)=9.8d0*rho_avg_global(tt)*head_pump(tt)       !averaged
            dp_err(tt)=dp_loop_global(tt)+dp_pump(tt)                !averaged 
            flux_err=cf*dp_err(tt)                                   !m3/s based on one face area or a total area by what the definition of cf is
            vol_flow(tt)=vol_flow(tt)+relax_flow*flux_err        !vol_flow and flux_err should be based on the same concept (one face or total face)
         ENDIF
!         
      ENDDO
!    
      IF(initial_write) THEN
         IF(myrank.eq.0) OPEN(500,file='MCP_data.dat')
         icnt_time=1
         write_interval=0.2d0  !user input: writing interval
         initial_write=.false.
      ENDIF
      IF(time.gt.write_interval*icnt_time) THEN 
         IF(myrank.eq.0) WRITE(500,102) time,(dp_pump(tt),tt=1,num_mcp),(vol_flow(tt),tt=1,num_mcp), &
                     (rho_avg_global(tt),tt=1,num_mcp),(head_mult(tt),tt=1,num_mcp),&
                     (head_pump(tt),tt=1,num_mcp),(dp_err(tt),tt=1,num_mcp), &
                     (rpm_controlled(tt),tt=1,num_mcp)
         icnt_time=icnt_time+1 
      ENDIF   
  102 FORMAT(100(e14.7,1x))
!    
!            
    END SUBROUTINE mcp_model
!    
!
!    
    SUBROUTINE mcp_head(hratio,spdrat_in,wratio_in)
!
!.....brief Read MCP head from the tabulated homologous curves.
!
!.....Returns normalized head, (head)/(rated head).
!
      USE Zmcp           ,ONLY: frac_tabl,han,hvn,had,hvd,hat,hvt,har,hvr  
!
      IMPLICIT NONE
!
      REAL(8),INTENT(in)::spdrat_in    !< normalized speed.
      REAL(8),INTENT(in)::wratio_in    !< normalized volume flow.
      REAL(8),INTENT(out)::hratio      !< normalized head from the homologous curves.
!
      REAL(8),PARAMETER::tol=1.d-9    !< small value threshold.
      REAL(8)::wratio                   !< (vol. flow)/(rated vol. flow)
      REAL(8)::spdrat                   !< (speed)/(rated speed)
      REAL(8)::pr                       !< a/v or v/a for selected  curve.
      REAL(8)::v2,a2                    !< squared volume flow ratio and speed ratio.
      INTEGER::regime                   !< homologous for current running condition.
      INTEGER::lpospr=1                 !< index in a homologous table.
!
      IF(DABS(spdrat_in).le.tol)THEN
         spdrat=0.0d0
      ELSE
         spdrat=spdrat_in
      ENDIF
!
      IF(DABS(wratio_in).le.tol)THEN
         wratio=0.0d0
      ELSE
         wratio=wratio_in
      ENDIF
!
      a2=spdrat*spdrat
      v2=wratio*wratio
!
!.....Search homologous curve to use.
!
      IF(spdrat.gt.0.0d0)THEN
         IF (wratio.gt.0.0d0)THEN
            IF(wratio.le.spdrat)THEN
               regime=1     !han
            ELSE
               regime=2     !hvn
            ENDIF
         ELSEIF(wratio.lt.0.0d0)THEN
            IF (wratio.ge.-spdrat)THEN
               regime=3     !had
            ELSE
               regime=4     !hvd
            ENDIF
         ELSE
            regime=1        !han
         ENDIF
      ELSEIF(spdrat.lt.0.0d0)THEN
         IF(wratio.lt.0.0d0)THEN
            IF(wratio.ge.spdrat) THEN
               regime=5     !hat
            ELSE
               regime=6     !hvt
            ENDIF
         ELSEIF(wratio.gt.0.0d0)THEN
            IF (wratio.le.-spdrat)THEN
               regime=7     !har
            ELSE
               regime=8     !hvr
            ENDIF
         ELSE
            regime=7        !har
         ENDIF
      ELSE
         IF(wratio.gt.0.0d0) THEN
            regime=2        !hvn
         ELSEIF(wratio.lt.0.0d0) THEN
            regime=4        !hvd
         ELSE
            regime=0
         ENDIF
      ENDIF
!
!.....Interpolate the selected curve.
!
      SELECT CASE(regime)
         CASE(1)  !han
            pr=wratio/spdrat
            CALL mcp_interp1(11,pr,frac_tabl,han,lpospr,hratio)
            hratio=hratio*a2
         CASE(2)  !hvn
            pr=spdrat/wratio
            CALL mcp_interp1(11,pr,frac_tabl,hvn,lpospr,hratio)
            hratio=hratio*v2
         CASE(3)  !had
            pr=1.0+wratio/spdrat
            CALL mcp_interp1(11,pr,frac_tabl,had,lpospr,hratio)
            hratio=hratio*a2
         CASE(4)  !hvd
            pr=1.0+spdrat/wratio
            CALL mcp_interp1(11,pr,frac_tabl,hvd,lpospr,hratio)
            hratio=hratio*v2
         CASE(5)  !hat
            pr=wratio/spdrat
            CALL mcp_interp1(11,pr,frac_tabl,hat,lpospr,hratio)
            hratio=hratio*a2
         CASE(6)  !hvt
            pr=spdrat/wratio
            CALL mcp_interp1(11,pr,frac_tabl,hvt,lpospr,hratio)
            hratio=hratio*v2
         CASE(7)  !har
            pr=1.0+wratio/spdrat
            CALL mcp_interp1(11,pr,frac_tabl,har,lpospr,hratio)
            hratio=hratio*a2
         CASE(8)  !hvr
            pr=1.0+spdrat/wratio
            CALL mcp_interp1(11,pr,frac_tabl,hvr,lpospr,hratio)
            hratio=hratio*v2
         CASE DEFAULT
            hratio=0.0d0
      END SELECT
!
      RETURN
      END SUBROUTINE mcp_head

!  
!    
      SUBROUTINE mcp_interp1(N, X, XTBL, YTBL, IPOS0, Y)
      !
      ! Single table linear interpolation routine.
      ! No extrapolation: If X is out of range, the table's end value will be returned.
      !
      ! I/O parameter
      !  n        length of array
      !  x        interpolation location
      !  xtbl     x-table
      !  ytbl     y-table
      !  ipos0    initial guess for location
      !           if 0, table is equi-spac, if non-zero, returns the actual location
      !  y        interpolated value
      !
      ! Affected external variables
      !  None
      !
      ! Local variables
      !  fact, ipos  dummy
      !
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: N
      REAL, INTENT(IN) :: X
      REAL, INTENT(IN), DIMENSION(N) :: XTBL, YTBL
      INTEGER, INTENT(INOUT) :: IPOS0
      REAL, INTENT(OUT) :: Y

      REAL :: FACT
      INTEGER :: IPOS

      IF ( N > 0 ) THEN
         CALL mcp_LOOKUP (N, X, XTBL, IPOS0, IPOS)
         IF (IPOS > 0) THEN
            FACT = (X - XTBL (IPOS) ) / (XTBL (IPOS + 1) - XTBL (IPOS) )
            Y    = YTBL (IPOS) + FACT * (YTBL (IPOS + 1) - YTBL (IPOS) )
         ELSE
            ! X is out of range.
            Y = YTBL (IABS (IPOS) )
         ENDIF
      else
         Y = -huge(Y)
      END IF

      RETURN

    END SUBROUTINE mcp_interp1    
    
    
    SUBROUTINE mcp_LOOKUP (N, X, XTBL, IPOS0, LOC)
      !
      ! Table look-up. Locates value X in the table XTBL.
      !
      ! I/O parameter
      !  n        length of array
      !  x        interpolation location
      !  xtbl     table
      !  ipos0    initial guess for location
      !           if 0, table is equi-spac, if non-zero, returns the actual location
      !  loc      located position. Returns negative if X is out of range
      !
      ! Affected external variables
      !  None
      !
      ! Local variables
      !  None
      !
      IMPLICIT NONE
   
      INTEGER, INTENT(IN) :: N
      REAL, INTENT(IN) :: X
      REAL, INTENT(IN), DIMENSION(N) :: XTBL
      INTEGER, INTENT(INOUT) :: IPOS0
      INTEGER, INTENT(OUT) :: LOC
   
      LOC = MIN (IPOS0, N - 1)
      IF (N == 1) THEN
         LOC = -1
         IPOS0 = 1
      ELSE IF ( (XTBL (N) - X) * (X - XTBL (1) ) < 0.0 ) THEN
         ! Check for X being out of range, to prevent extrapolation.
         LOC = - N
         IF ( (XTBL (2) - XTBL (1) ) * (XTBL (1) - X) > 0.0 ) then
            LOC = - 1
         end if
         IF ( IPOS0 > 0 ) then
            IPOS0 = - LOC
         end if
   
      ELSEIF (IPOS0 <= 0) THEN
         ! Calculate position in equi-spaced XTBL:
         LOC = IFIX ( (X - XTBL (1) ) / (XTBL (2) - XTBL (1) ) ) + 1
         LOC = MIN  (N - 1, LOC)
   
      ELSEIF (XTBL (2) > XTBL (1) ) THEN
         ! Find position in vari-spaced XTBL, monotonically increasing.
         IF (X < XTBL (LOC) ) THEN
            LOC = LOC - 1
            DO WHILE ( X < XTBL (LOC) )
               LOC = LOC - 1
            END DO
         ELSE
            DO WHILE ( X > XTBL (LOC + 1) )
               LOC = LOC + 1
            END DO
         ENDIF
         IPOS0 = LOC
   
      ELSE
         ! Find position in vari-spaced XTBL, monotonically decreasing.
         IF (X > XTBL (LOC) ) THEN
            LOC = LOC - 1
            DO WHILE (X > XTBL (LOC) )
               LOC = LOC - 1
            END DO
         ELSE
            DO WHILE (X < XTBL (LOC + 1) )
               LOC = LOC + 1
            END DO
         ENDIF
         IPOS0 = LOC
      ENDIF
      END SUBROUTINE mcp_LOOKUP
