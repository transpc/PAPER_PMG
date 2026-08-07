!   
      SUBROUTINE radiation_component
!      
      USE Zzone      ,ONLY: ncell_fluid
      USE Vol_DATA   ,ONLY: cell
      USE Zwall_HTC  ,ONLY: twall_rv
      USE Solid_DATA ,ONLY: solid
      USE Zrad_comp 
      USE Zcore      ,ONLY: np
!
      USE Zrv_ncell  ,ONLY: ncell_fluid_core,cupid_cell_channel,n_channel_fluid,nz_fluid
!      
      IMPLICIT NONE
!      
      INTEGER i,j,k,m
      LOGICAL,SAVE::initial=.TRUE.
      REAL(8)::htc
      REAL(8),ALLOCATABLE:: tary(:),tary_tmp(:)
!      
!.....Start to calculate radiation 
!      
      IF(input_opt.eq.3)THEN
         CALL radiation_component_general
         RETURN
      ENDIF   
!      
      IF(initial)THEN
!         
         initial=.FALSE.
         ALLOCATE(heat_rad(nsize_max,nset))
         ALLOCATE(t1(nsize_max,nset))
         ALLOCATE(t2(nsize_max,nset))
         ALLOCATE(heat_rad1(nsize_max,nset))
         ALLOCATE(heat_rad2(nsize_max,nset))
         ALLOCATE(rho(2,nset))
         ALLOCATE(a(nset),b(nset),c(nset),d(nset),dterm(nset))
!         
         heat_rad=0.0d0
         t1=0.0d0
         t2=0.0d0
         heat_rad1=0.0d0
         heat_rad2=0.0d0
         rho=0.0d0
         a=0.0d0
         b=0.0d0
         c=0.0d0
         d=0.0d0
         dterm=0.0d0
!
!.....ncell_fluid_core index for twall_rv(:,1)
         !ALLOCATE(core_cell_cupid(ncell_fluid_core))
         ALLOCATE(core_cell_cupid(ncell_fluid))  
         core_cell_cupid=0
         DO i=1,ncell_fluid_core
            j=n_channel_fluid(i)
            k=nz_fluid(i)
            m=cupid_cell_channel(i) !m=cupid_cell, i=core cell
            core_cell_cupid(m)=i    
         ENDDO
         DO j=1,nset
            IF(nproperty(1,j).eq.1)THEN
               DO i=1,nsize(1,j)
                   IF(cell_idx(i,1,j).eq.0)CYCLE 
                   k=core_cell_cupid(cell_idx(i,1,j)) 
                   IF(k.ge.1.and.k.le.ncell_fluid_core)THEN
                      core_idx(i,1,j)=k
                   ELSE
                      WRITE(*,*)'1core_idx error!'
                      PAUSE
                      STOP
                   ENDIF   
               ENDDO
            ENDIF
            IF(nproperty(2,j).eq.1)THEN
               DO i=1,nsize(1,j)   
                   IF(cell_idx(i,2,j).eq.0)CYCLE 
                   k=core_cell_cupid(cell_idx(i,2,j)) 
                   IF(k.ge.1.and.k.le.ncell_fluid_core)THEN
                      core_idx(i,2,j)=k             
                   ELSE
                      WRITE(*,*)'2core_idx error!'
                      PAUSE
                      STOP
                   ENDIF   
               ENDDO
            ENDIF   
         ENDDO   
         DEALLOCATE(core_cell_cupid)
!
!........radiation heat transfer maxtrix coefficients
      IF(1)THEN         
         DO j=1,nset
            rho(:,j)=1.0d0-epsil(:,j)
            a(j)=1.0d0-rho(1,j)*viewf(1,1,j)
            b(j)=-rho(1,j)*viewf(1,2,j)
            c(j)=-rho(2,j)*viewf(2,1,j)
            d(j)=1.0d0-rho(2,j)*viewf(2,2,j)            
            dterm(j)=a(j)*d(j)-b(j)*c(j)
         ENDDO
      ELSE         
         DO j=1,nset
            rho(:,j)=1.0d0-epsil(:,j)
            a(j)=1.0d0-rho(1,j)*viewf(1,1,j)
            b(j)=-rho(1,j)*viewf(2,1,j)
            c(j)=-rho(2,j)*viewf(1,2,j)
            d(j)=1.0d0-rho(2,j)*viewf(2,2,j)            
            dterm(j)=a(j)*d(j)-b(j)*c(j)
         ENDDO         
      ENDIF         
!         
         RETURN
      ENDIF
!      
      heat_rad=0.0d0
      qrad_flu=0.0d0
      qrad_rod=0.0d0
      qrad_sol=0.0d0
      t1=0.0d0
      t2=0.0d0
      heat_rad1=0.0d0
      heat_rad2=0.0d0      
!
!.....Calculate radiation heat flux
!    
!
!.....pick up t1, t2 at each core      
      DO j=1,nset
         DO i=1,nsize(1,j)
            k=cell_idx(i,1,j)
            IF(k.eq.0)CYCLE
            IF(nproperty(1,j).eq.0)THEN !fluid
               t1(i,j)=cell%tl(k)
            ELSEIF(nproperty(1,j).eq.1)THEN !fuel rod in fluid
               k=core_idx(i,1,j)
               t1(i,j)=twall_rv(k,1) !nr_2d
            ELSEIF(nproperty(1,j).eq.2)THEN !pure solid or solid in porous media
               t1(i,j)=solid%tsol(k)
            ELSE
               WRITE(*,*)'error1 in radiation_component!'
               PAUSE
               STOP
            ENDIF   
         ENDDO  
      ENDDO              
!            
      DO j=1,nset
         DO i=1,nsize(1,j)
            k=cell_idx(i,2,j)
            IF(k.eq.0)CYCLE
            IF(nproperty(2,j).eq.0)THEN
               t2(i,j)=cell%tl(k)
            ELSEIF(nproperty(2,j).eq.1)THEN
               k=core_idx(i,2,j)
               t2(i,j)=twall_rv(k,1) 
            ELSEIF(nproperty(2,j).eq.2)THEN
               t2(i,j)=solid%tsol(k)
            ELSE   
               WRITE(*,*)'error2 in radiation_component!'
               PAUSE
               STOP               
            ENDIF             
         ENDDO  
      ENDDO          
!
!.....complete t1, t2 for all the cores      
      IF(np.gt.1)THEN
         ALLOCATE(tary(nsize_max),tary_tmp(nsize_max))
         DO j=1,nset
            tary(:)=t1(:,j)
            CALL allreduce_r(tary,tary_tmp,nsize(1,j))
            t1(:,j)=tary_tmp(:)
            tary(:)=t2(:,j)
            CALL allreduce_r(tary,tary_tmp,nsize(1,j))
            t2(:,j)=tary_tmp(:)            
         ENDDO 
         DEALLOCATE(tary,tary_tmp)
      ENDIF
!   
!.....Calculate radiation heat flux for all the cores     
      DO j=1,nset
         DO i=1,nsize(1,j)      
            htc=sigma_sb*epsil(1,j)*epsil(2,j)*area(1,j)*(t1(i,j)+t2(i,j))*(t1(i,j)*t1(i,j)+t2(i,j)*t2(i,j))
            heat_rad(i,j)=htc*(t1(i,j)-t2(i,j)) 
            heat_rad1(i,j)=+d(j)/dterm(j)*epsil(1,j)*sigma_sb*t1(i,j)**4.0d0-b(j)/dterm(j)*epsil(2,j)*sigma_sb*t2(i,j)**4.0d0
            heat_rad2(i,j)=-c(j)/dterm(j)*epsil(1,j)*sigma_sb*t1(i,j)**4.0d0+a(j)/dterm(j)*epsil(2,j)*sigma_sb*t2(i,j)**4.0d0
            heat_rad1(i,j)=epsil(1,j)/rho(1,j)*(sigma_sb*t1(i,j)**4.0d0-heat_rad1(i,j))
            heat_rad2(i,j)=epsil(2,j)/rho(2,j)*(sigma_sb*t2(i,j)**4.0d0-heat_rad2(i,j))
            heat_rad1(i,j)=heat_rad1(i,j)*area(1,j)
            heat_rad2(i,j)=heat_rad2(i,j)*area(2,j)
            heat_rad(i,j)=heat_rad1(i,j)
         ENDDO  
      ENDDO
!
!.....Distribute radiation heat flux to the source term of fluid, rod, solid to each core
!      
      IF(0)THEN      
         DO j=1,nset
            DO i=1,nsize(1,j)
               k=cell_idx(i,1,j)
               IF(k.eq.0)CYCLE
               IF(nproperty(1,j).eq.0)THEN !fluid
                  qrad_flu(k)=qrad_flu(k)-heat_rad(i,j)
               ELSEIF(nproperty(1,j).eq.1)THEN !fuel rod in fluid
                  qrad_rod(k)=qrad_rod(k)-heat_rad(i,j)
               ELSEIF(nproperty(1,j).eq.2)THEN !pure solid or solid in porous media
                  qrad_sol(k)=qrad_sol(k)-heat_rad(i,j)
               ENDIF   
            ENDDO  
         ENDDO               
         DO j=1,nset
            DO i=1,nsize(1,j)                
               k=cell_idx(i,2,j)
               IF(k.eq.0)CYCLE
               IF(nproperty(2,j).eq.0)THEN
                  qrad_flu(k)=qrad_flu(k)+heat_rad(i,j)
               ELSEIF(nproperty(2,j).eq.1)THEN
                  qrad_rod(k)=qrad_rod(k)+heat_rad(i,j)
               ELSEIF(nproperty(2,j).eq.2)THEN
                  qrad_sol(k)= qrad_sol(k)+heat_rad(i,j)
               ENDIF             
            ENDDO  
         ENDDO  
      ELSE  
         DO j=1,nset
            DO i=1,nsize(1,j)
               k=cell_idx(i,1,j)
               IF(k.eq.0)CYCLE
               IF(nproperty(1,j).eq.0)THEN !fluid
                  qrad_flu(k)=qrad_flu(k)-heat_rad1(i,j)
               ELSEIF(nproperty(1,j).eq.1)THEN !fuel rod in fluid
                  qrad_rod(k)=qrad_rod(k)-heat_rad1(i,j)
               ELSEIF(nproperty(1,j).eq.2)THEN !pure solid or solid in porous media
                  qrad_sol(k)=qrad_sol(k)-heat_rad1(i,j)
               ENDIF   
            ENDDO  
         ENDDO               
         DO j=1,nset
            DO i=1,nsize(1,j)         
               k=cell_idx(i,2,j)
               IF(k.eq.0)CYCLE
               IF(nproperty(2,j).eq.0)THEN
                  qrad_flu(k)=qrad_flu(k)-heat_rad2(i,j)
               ELSEIF(nproperty(2,j).eq.1)THEN
                  qrad_rod(k)=qrad_rod(k)-heat_rad2(i,j)
               ELSEIF(nproperty(2,j).eq.2)THEN
                  qrad_sol(k)= qrad_sol(k)-heat_rad2(i,j)
               ENDIF             
            ENDDO  
         ENDDO
      ENDIF      
!
      RETURN
      ENDSUBROUTINE radiation_component
!   
