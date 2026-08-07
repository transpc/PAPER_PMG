!   
      SUBROUTINE radiation_component_init(ns_n,nn,nzone_tmp,xloc_tmp) !pik-radiation_component
!
      USE Zparam     ,ONLY: ndim
      USE Zrad_comp  
!      
      IMPLICIT NONE
!      
      INTEGER:: nn,ns_n(nn) !the geometrical number of cells , solid cell index
      INTEGER:: nzone_tmp(nn)
      REAL(8):: xloc_tmp(nn,ndim)
!       
      INTEGER:: i,j
      INTEGER:: err
!      
      OPEN(333,file='radiation_component.in',status='old',iostat=err)
      IF(err.ne.0)THEN
         rad_comp_mod=0
         RETURN        
      ELSE
         rad_comp_mod=1
      ENDIF   
!      
      READ(333,*)input_opt
!      
      IF(input_opt.eq.1)THEN
         CLOSE(333)
         CALL radiation_component_init_zone(nn,nzone_tmp,xloc_tmp)
      ELSEIF(input_opt.eq.3)THEN
         CLOSE(333)
         CALL radiation_component_init_general
      ELSE      
!
!.....Read radiation information
!     
         READ(333,*)nset,nsize_max
         ALLOCATE(nproperty(2,nset))
         ALLOCATE(nsize(2,nset))
         ALLOCATE(epsil(2,nset))
         ALLOCATE(viewf(2,2,nset))      
         ALLOCATE(area(2,nset))
         ALLOCATE(cell_idx_tmp(nsize_max,2,nset))
         nproperty=0
         nsize=0
         epsil=0.0d0
         viewf=0.0d0
         area=0.0d0
         cell_idx_tmp(:,:,:)=0
!      
         DO j=1,nset
            READ(333,*)nproperty(1,j),nproperty(2,j),nsize(1,j)              !fluid solid size 
            nsize(2,j)=nsize(1,j)
            READ(333,*)epsil(1,j),epsil(2,j),area(1,j),area(2,j)           !epilson, radiation, area
            READ(333,*)viewf(1,1,j),viewf(1,2,j),viewf(2,1,j),viewf(2,2,j) !view factor
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
      ENDIF !input_opt  
!
      CALL radiation_component_init_parallel(ns_n,nn)
!
      END SUBROUTINE radiation_component_init
!   
!------------------------------------------------------------------------------      
!------------------------------------------------------------------------------      
!      
      SUBROUTINE radiation_component_init_parallel(ns_n,nn)
!
      USE Zzone      ,ONLY: ncell_fluid,ncell_cond
      USE Zmpi       ,ONLY: jperm,jperm_c
      USE Zcore      ,ONLY: myrank,np
      USE Zrad_comp  
!      
      IMPLICIT NONE
!      
      INTEGER:: nn,ns_n(nn) !the geometrical number of cells , solid cell index
      INTEGER:: i,j,k
      INTEGER:: icell,check
      INTEGER:: idx1,idx2      

!
!.....For now, we have common part of radiation_component_init and radiation_component_init_cell
!
!
!.....convert raw cell index into 1~ncell_fluid_all and 1~ncell_cond_all
!      
      IF(myrank.eq.0)THEN         
         DO j=1,nset
            DO i=1,nsize(1,j)
               idx1=cell_idx_tmp(i,1,j)
               idx2=cell_idx_tmp(i,2,j)
               IF(nproperty(1,j).eq.0)THEN !fluid
                  cell_idx_tmp(i,1,j)=idx1
               ELSEIF(nproperty(1,j).eq.1)THEN !fuel rod in fluid
                  cell_idx_tmp(i,1,j)=idx1
               ELSEIF(nproperty(1,j).eq.2)THEN !pure solid or solid in porous media
                  cell_idx_tmp(i,1,j)=ns_n(idx1)
               ENDIF   
               IF(nproperty(2,j).eq.0)THEN !fluid
                  cell_idx_tmp(i,2,j)=idx2
               ELSEIF(nproperty(2,j).eq.1)THEN !fuel rod in fluid
                  cell_idx_tmp(i,2,j)=idx2
               ELSEIF(nproperty(2,j).eq.2)THEN !pure solid or solid in porous media
                  cell_idx_tmp(i,2,j)=ns_n(idx2)
               ENDIF
            ENDDO
         ENDDO 
!
!........For special output         
         DO i=1,ncell_out
            idx1=cell_out_tmp(i,1)
            idx2=cell_out_tmp(i,2)
            IF(nproperty_out_tmp(i).eq.0)THEN !fluid
               cell_out_tmp(i,1)=idx1
            ELSEIF(nproperty_out_tmp(i).eq.1)THEN !fuel rod in fluid
               cell_out_tmp(i,1)=idx1
            ELSEIF(nproperty_out_tmp(i).eq.2)THEN !pure solid or solid in porous media
               cell_out_tmp(i,1)=ns_n(idx1)
            ENDIF 
            IF(nproperty_out_tmp(i).eq.0)THEN !fluid
               cell_out_tmp(i,2)=idx2
            ELSEIF(nproperty_out_tmp(i).eq.1)THEN !fuel rod in fluid
               cell_out_tmp(i,2)=idx2
            ELSEIF(nproperty_out_tmp(i).eq.2)THEN !pure solid or solid in porous media
               cell_out_tmp(i,2)=ns_n(idx2)
            ENDIF             
         ENDDO   
!         
      ENDIF !myrank.eq.0 
!
!      
      IF(np.gt.1) THEN
         DO j=1,nset
            CALL broadcast_i(cell_idx_tmp(:,1,j),nsize_max)
            CALL broadcast_i(cell_idx_tmp(:,2,j),nsize_max)
         ENDDO  
         CALL broadcast_i(cell_out_tmp(:,1),ncell_out)
         CALL broadcast_i(cell_out_tmp(:,2),ncell_out)
         CALL broadcast_i(nproperty_out_tmp(:),ncell_out)         
         CALL broadcast_r(frac_out_tmp(:,1),ncell_out)
         CALL broadcast_r(frac_out_tmp(:,2),ncell_out)
      ENDIF
!
!.....domain decomposition of cell_idx_tmp
!      
      ALLOCATE(cell_idx(nsize_max,2,nset))
      ALLOCATE(core_idx(nsize_max,2,nset))
      cell_idx=0
      core_idx=0
      DO j=1,nset
         DO i=1,nsize(1,j)
            check=0
            DO icell=1,ncell_fluid
               k=cell_idx_tmp(i,1,j)
               IF(nproperty(1,j).eq.0.or.nproperty(1,j).eq.1)THEN 
                  IF(jperm(icell).eq.k)THEN
                     cell_idx(i,1,j)=icell
                     check=check+1
                  ENDIF   
               ENDIF   
               k=cell_idx_tmp(i,2,j)
               IF(nproperty(2,j).eq.0.or.nproperty(2,j).eq.1)THEN
                  IF(jperm(icell).eq.k)THEN
                     cell_idx(i,2,j)=icell
                     check=check+1
                  ENDIF   
               ENDIF             
            ENDDO !icell
!
            DO icell=1,ncell_cond
               k=cell_idx_tmp(i,1,j) 
               IF(nproperty(1,j).eq.2)THEN 
                  IF(jperm_c(icell).eq.k)THEN
                     cell_idx(i,1,j)=icell !is this right?
                     check=check+1
                  ENDIF   
               ENDIF   
               k=cell_idx_tmp(i,2,j) 
               IF(nproperty(2,j).eq.2)THEN
                  IF(jperm_c(icell).eq.k)THEN
                     cell_idx(i,2,j)=icell
                     check=check+1
                  ENDIF   
               ENDIF             
            ENDDO !icell   
            IF(check.ne.2.and.check.ne.0)THEN
               !WRITE(*,"(a,1i3)")'No pair in radiation_component.in!!!,check=',check
               !PAUSE
               !STOP
            ENDIF   
         ENDDO !i 
      ENDDO !j 
!
!.....Print results
!      
      IF(0)THEN      
         DO j=1,nset
            WRITE(*,"(a,20i3)")'cell_idx1=',myrank,(cell_idx(i,1,nset),i=1,nsize_max),nproperty(1,j)
            WRITE(*,"(a,20i3)")'cell_idx2=',myrank,(cell_idx(i,2,nset),i=1,nsize_max),nproperty(2,j)
         ENDDO
         WRITE(*,"(a,20i3)")'cell_out1=',(cell_out_tmp(i,1),i=1,ncell_out)
         WRITE(*,"(a,20i3)")'cell_out2=',(cell_out_tmp(i,2),i=1,ncell_out)
         WRITE(*,"(a,20i3)")'nproperty=',(nproperty_out_tmp(i),i=1,ncell_out)
         WRITE(*,"(a,20f8.3)")'frac_out=',(frac_out_tmp(i,1),i=1,ncell_out)
         WRITE(*,"(a,20f8.3)")'frac_out=',(frac_out_tmp(i,2),i=1,ncell_out)
      ENDIF      
!
      END SUBROUTINE radiation_component_init_parallel
!   
