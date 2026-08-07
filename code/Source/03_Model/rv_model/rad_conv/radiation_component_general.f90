!   
      SUBROUTINE radiation_component_init_general !pik-radiation_component
!
      USE Zrad_comp  
!      
      IMPLICIT NONE
!       
      INTEGER:: i,j,k
!
!.....Read radiation information
!     
      OPEN(333,file='radiation_component.in',status='old')
      READ(333,*)input_opt
      READ(333,*)nset,nsize_max
      ALLOCATE(nproperty(2,nset))
      ALLOCATE(nsize(2,nset))
      ALLOCATE(epsil(nsize_max*2,nset))
      ALLOCATE(area(nsize_max*2,nset))
      ALLOCATE(viewf(nsize_max*2,nsize_max*2,nset))      
      ALLOCATE(cell_idx_tmp(nsize_max,2,nset))
      nproperty=0
      nsize=0
      epsil=0.0d0
      area=0.0d0
      cell_idx_tmp=0
      DO j=1,nset
         READ(333,*)nproperty(1,j),nproperty(2,j),nsize(1,j),nsize(2,j)              !fluid solid size 
         READ(333,*)(epsil(k,j),k=1,nsize(1,j))
         READ(333,*)(epsil(k+nsize(1,j),j),k=1,nsize(2,j))
         READ(333,*)(area(k,j),k=1,nsize(1,j))
         READ(333,*)(area(k+nsize(1,j),j),k=1,nsize(2,j))         
         DO i=1,nsize(1,j)  !mine-1,1(2사)
            READ(333,*)(viewf(i,k,j),k=1,nsize(1,j))
         ENDDO    
         DO i=1,nsize(1,j)  !partner-1,2(1사)
            READ(333,*)(viewf(i,nsize(1,j)+k,j),k=1,nsize(2,j))
         ENDDO   
         DO i=1,nsize(2,j)  !mine-2,2(4사)
            READ(333,*)(viewf(nsize(1,j)+i,nsize(1,j)+k,j),k=1,nsize(2,j))
         ENDDO 
         DO i=1,nsize(2,j)  !partner-2,1(3사)
            READ(333,*)(viewf(nsize(1,j)+i,k,j),k=1,nsize(1,j))
         ENDDO          
         READ(333,*)(cell_idx_tmp(i,1,j),i=1,nsize(1,j))
         READ(333,*)(cell_idx_tmp(i,2,j),i=1,nsize(2,j))
      ENDDO
!
!.....For special output      
      READ(333,*)ncell_out
      ALLOCATE(cell_out_tmp(ncell_out,2))
      ALLOCATE(frac_out_tmp(ncell_out,2))
      ALLOCATE(nproperty_out_tmp(ncell_out))
      DO i=1,ncell_out
         READ(333,*)cell_out_tmp(i,1),cell_out_tmp(i,2),frac_out_tmp(i,1),frac_out_tmp(i,2),nproperty_out_tmp(i)
      ENDDO   
!      
      CLOSE(333)
!
      END SUBROUTINE radiation_component_init_general
!   
!------------------------------------------------------------------------------      
!------------------------------------------------------------------------------      
!   
      SUBROUTINE radiation_component_general
!      
      USE Zzone      ,ONLY: ncell_fluid
      USE Vol_DATA   ,ONLY: cell
      USE Zwall_HTC  ,ONLY: twall_rv
      USE Solid_DATA ,ONLY: solid
      USE Zrad_comp 
!
      USE Zrv_ncell  ,ONLY: ncell_fluid_core,cupid_cell_channel,n_channel_fluid,nz_fluid
!      
      IMPLICIT NONE
!      
      INTEGER i,j,k,m
      LOGICAL,SAVE::initial=.TRUE.
      REAL(8),ALLOCATABLE,SAVE::am(:,:,:),bm(:,:)
      INTEGER,ALLOCATABLE,SAVE::pivot(:,:)
      REAL(8),ALLOCATABLE::delta(:,:)
      !REAL(8)::t1,t2
!      
      IF(initial)THEN
!         
         initial=.FALSE.
         ALLOCATE(heat_rad(nsize_max,nset))
         ALLOCATE(t1(nsize_max,nset))
         ALLOCATE(t2(nsize_max,nset))
         ALLOCATE(heat_rad1(nsize_max,nset))
         ALLOCATE(heat_rad2(nsize_max,nset))
         ALLOCATE(rho(nsize_max*2,nset))
         ALLOCATE(am(nsize_max*2,nsize_max*2,nset))
         ALLOCATE(bm(nsize_max*2,nset))
         ALLOCATE(pivot(nsize_max*2,nset))
         DO j=1,nset
            rho(:,j)=1.0d0-epsil(:,j)
         ENDDO   
         heat_rad=0.0d0
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
                      WRITE(*,*)'3core_idx error!'
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
                      WRITE(*,*)'4core_idx error!'
                      PAUSE
                      STOP
                   ENDIF   
               ENDDO
            ENDIF   
         ENDDO   
         DEALLOCATE(core_cell_cupid)
!
!........radiation heat transfer maxtrix coefficients
         ALLOCATE(delta(nsize_max*2,nsize_max*2))
         delta(:,:)=0.0d0
         DO i=1,nsize_max*2
           delta(i,i)=1.0d0
         ENDDO  
         DO j=1,nset
            DO i=1,nsize(1,j)+nsize(2,j)
               DO k=1,nsize(1,j)+nsize(2,j)
                  am(i,k,j)=delta(i,k)-rho(i,j)*viewf(i,k,j)
               ENDDO   
            ENDDO   
         ENDDO
         DEALLOCATE(delta)
!
         pivot(:,:)=0
         DO j=1,nset
            call ddgef_rad(am(:,:,j),nsize(1,j)+nsize(2,j),pivot(:,j))
         ENDDO   
!         
         RETURN
      ENDIF
!      
      heat_rad=0.0d0
      qrad_flu=0.0d0
      qrad_rod=0.0d0
      qrad_sol=0.0d0
!
!.....Calculate radiation heat flux
!      
      DO j=1,nset
         DO i=1,nsize(1,j)
            k=cell_idx(i,1,j)
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
            bm(i,j)=epsil(i,j)*sigma_sb*t1(i,j)**4.0d0
         ENDDO  
         DO i=1,nsize(2,j)
            k=cell_idx(i,2,j)
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
            bm(i+nsize(1,j),j)=epsil(i+nsize(1,j),j)*sigma_sb*t2(i,j)**4.0d0
         ENDDO  
!         
         CALL ddges_rad(am(:,:,j),nsize(1,j)+nsize(2,j),pivot(:,j),bm(:,j))            
!         
         DO i=1,nsize(1,j)
            heat_rad1(i,j)=bm(i,j)
         ENDDO   
         DO i=1,nsize(2,j)
            k=i+nsize(1,j)            
            heat_rad2(i,j)=bm(k,j)
         ENDDO
!        
         DO i=1,nsize(1,j)
            heat_rad1(i,j)=epsil(i,j)/rho(i,j)*(sigma_sb*t1(i,j)**4.0d0-heat_rad1(i,j))
            heat_rad1(i,j)=heat_rad1(i,j)*area(i,j)
         ENDDO
         DO i=1,nsize(2,j)
            k=i+nsize(1,j)
            heat_rad2(i,j)=epsil(k,j)/rho(k,j)*(sigma_sb*t2(i,j)**4.0d0-heat_rad2(i,j))
            heat_rad2(i,j)=heat_rad2(i,j)*area(k,j)         
         ENDDO   
      ENDDO
!
!.....Distribute radiation heat flux to the source term of fluid, rod, solid
!      
      IF(0)THEN      
         DO j=1,nset
            DO i=1,nsize(1,j)
               k=cell_idx(i,1,j)
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
      END SUBROUTINE radiation_component_general
!----------------------------------------------------------------------
      SUBROUTINE ddgef_rad(a,n,pivot)
!
      IMPLICIT NONE
!.....Input
      INTEGER :: n
      INTEGER,DIMENSION(n) :: pivot
      REAL(8),DIMENSION(n,n) :: a
!
      INTEGER :: i,j,k,ii,jj,ip
      REAL(8) :: temp,pivotx
!
       DO j=1,n
!-----find the pivot and pivot column j
          pivotx=a(j,j)
          ip=j
          do ii=j+1,n
             if(abs(pivotx).lt.abs(a(ii,j))) then
                pivotx=a(ii,j)
                ip=ii
             endif
          enddo
          pivot(j)=ip
          a(ip,j)=a(j,j)
          pivotx=1.d0/pivotx
          a(j,j)=pivotx
          do jj=j+1,n
             temp=a(ip,jj)
             a(ip,jj)=a(j,jj)
             a(j,jj)=temp
          enddo
!--------
          do i=j+1,n
             a(i,j)=a(i,j)*pivotx
             do k=j+1,n
                temp=a(i,k)
                temp=temp-a(i,j)*a(j,k)
                a(i,k)=temp
             enddo
          enddo
!
       enddo
!-------back pivot
        do jj=1,n
           do i=jj+1,n 
             ip=pivot(i)
             temp=a(ip,jj)
             a(ip,jj)=a(i,jj)
             a(i,jj)=temp
           enddo
        enddo
!
      END SUBROUTINE ddgef_rad
!   
       SUBROUTINE ddges_rad(a,n,pivot,b)
!
      IMPLICIT NONE
!
!.....Input
      INTEGER :: n
      INTEGER,DIMENSION(n) :: pivot
      REAL(8),DIMENSION(n,n) :: a
      REAL(8),DIMENSION(n) :: b
!
      INTEGER :: i,j,k,ip
      REAL(8) :: temp
!------
      do i=1,n
         ip=pivot(i) 
         temp=b(ip)
         b(ip)=b(i)
         b(i)=temp
      enddo

!------     
      do j=2,n
         temp=b(j)
         do k=1,j-1
            temp=temp-a(j,k)*b(k)
         enddo
         b(j)=temp
      enddo
!------     
      do j=n,1,-1
         temp=b(j)
         do k=j+1,n
            temp=temp-a(j,k)*b(k)
         enddo
         b(j)=temp*a(j,j)
       enddo
!
       END SUBROUTINE ddges_rad
   
