!
      SUBROUTINE valve_define_fluxBC
!
!     DATA STRUCTURE
!
!     num_valve=total number of valves in global domain
!     num_valveloc=number of valves in local subdomain
!     num_valveface(1:num_valveloc)=number of faces of a valve in local subdomain
!     n_face_valve(1:num_valveloc,1:num_valveface)=face index of a valve in local subdomain
!     fzone_valve(1:num_valve,1)=fluid zone index of valve in global domain (rv_parameters.in)
!     fzone_valve(1:num_valve,2)=fluid zone index of valve in global domain (rv_parameters.in)
!  
      USE Zinterface
      USE Zzone        , ONLY: nzone,ncell_fluid_all,ncell_fluid_all
      USE Znum_cell    , ONLY: istart_nf  
      USE Zvec_index   , ONLY: left_nf,right_non
      USE Zvec_param   , ONLY: nf_nonk
      USE Znum_cell    , ONLY: right_nb_k
      USE Zcore        , ONLY: np,myrank
      USE Zvalve       , ONLY: num_valveloc,num_valve,num_valveface,n_face_valve, &
                               mapping_valve,icell_valve,nonk_valve,fzone_valve,  &
                               ocell_valve
!
      IMPLICIT NONE
!
      INTEGER :: get_global_cell       
      INTEGER :: nf_number,istart,isize,i1,ir,j,k
      INTEGER :: i,ii,kk,tt,nn,zz,ii_valve,kk_valve,ii_tmp,kk_tmp
      INTEGER :: ivalve,ovalve
      INTEGER :: err,unit_valve,num_zone,iig,kkg,k1,k2
      INTEGER,DIMENSION(:),ALLOCATABLE :: cell_fluxBC,zone_fluxBC,nzone_fluxBC,nn_array
      LOGICAL,SAVE :: initial_loc
      LOGICAL,SAVE :: initial_nzone=.true.
!
!.....Open FluxBCzone file
!      
      unit_valve=422 
      OPEN(unit_valve,file='FluxBCZone.in',status='old',iostat=err)
      IF(err.eq.0)then
         IF(myrank.eq.0)WRITE(*       ,"(11x,a)")'Reading FluxBCZone_valve.in...'
         READ(unit_valve,*) num_zone
         ALLOCATE(cell_fluxBC(num_zone),zone_fluxBC(num_zone),nzone_fluxBC(ncell_fluid_all))
         nzone_fluxBC=0
         DO k=1,num_zone
            READ(unit_valve,*) cell_fluxBC(k),zone_fluxBC(k)
            ii=cell_fluxBC(k)
            nzone_fluxBC(ii)=zone_fluxBC(k) !ii=global, nzone_valve=FluxBC zone number
         ENDDO
      ELSE
         print*,'>> ERROR! FluxBCZone.in file must be in the folder.'
         pause
         stop
      ENDIF
      CLOSE(unit_valve)
!
!......Communicate nzone
!      
      IF(initial_nzone) THEN
         PRINT*,'>>>valve face is finding...'
         IF(np.gt.1) CALL communicate_1d_int(nzone)         
         initial_nzone=.false.
      ENDIF   
!
!......Initialize the array size
!       
      num_valveloc=0
      num_valveface=0
      nf_number=0
      istart=istart_nf(1,nf_number)
      isize =istart_nf(2,nf_number)
      DO tt=1,num_valve            !num_mcp=total number of mcps of a reactor        
         ivalve=fzone_valve(tt,1)     !inflow fluid zone number of a valve
         ovalve=fzone_valve(tt,2)     !outflow fluid zone number of a valve  
!
         initial_loc=.true. 
         DO i=1,isize
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)
            iig=get_global_cell(ii) !global
            kkg=get_global_cell(kk) !global            
!               
            k1=0
            k2=0
            DO k=1,num_zone
               IF(iig.eq.cell_fluxBC(k)) THEN
                  k1=k
               ENDIF
               IF(kkg.eq.cell_fluxBC(k)) THEN
                  k2=k
               ENDIF               
            ENDDO
            IF(k1.ne.0.and.k2.ne.0) THEN
               IF((zone_fluxBC(k1).eq.ivalve.and.zone_fluxBC(k2).eq.ovalve).or. &
                  (zone_fluxBC(k2).eq.ivalve.and.zone_fluxBC(k1).eq.ovalve)) THEN            
!                   
                  IF(initial_loc) THEN
                     num_valveloc=num_valveloc+1  
                     initial_loc=.false.
                  ENDIF  
                  IF(num_valveloc.ne.0) num_valveface(num_valveloc)=num_valveface(num_valveloc)+1               
               ENDIF
            ENDIF          
         ENDDO     
      ENDDO         
!
!......Initialize the array wrt the max value of valve faces
!      
      ALLOCATE(nn_array(np))
      nn_array=0
      nn_array(myrank+1)=MAXVAL(num_valveface(1:num_valve)) 
      CALL allreducei_i(nn_array,np)       
      nn=MAXVAL(nn_array)       
!      nn=MAXVAL(num_valveface(1:num_valve))
      IF(.not.ALLOCATED(n_face_valve)) ALLOCATE(n_face_valve(num_valve,nn))     !local domain 
      IF(.not.ALLOCATED(mapping_valve)) ALLOCATE(mapping_valve(num_valve))      !local-to-global
      IF(.not.ALLOCATED(icell_valve)) ALLOCATE(icell_valve(num_valve,nn))       !local domain
      IF(.not.ALLOCATED(ocell_valve)) ALLOCATE(ocell_valve(num_valve,nn))       !local domain
      IF(.not.ALLOCATED(nonk_valve)) ALLOCATE(nonk_valve(num_valve,nn))        !local domain       
!      
!......Finding valve faces
!       
      num_valveloc=0
      num_valveface=0
      nn=0
      istart=istart_nf(1,nf_number)
      isize =istart_nf(2,nf_number)          
      DO tt=1,num_valve 
         ivalve=fzone_valve(tt,1)     !inflow fluid zone number of a valve
         ovalve=fzone_valve(tt,2)     !outflow fluid zone number of a valve             
!             
         initial_loc=.true.
         DO i=1,isize
            i1=istart+i
            ii=left_nf(i1)   !local
            kk=right_non(i)  !local 
            iig=get_global_cell(ii) !global
            kkg=get_global_cell(kk) !global

            k1=0
            k2=0
            DO k=1,num_zone
               IF(iig.eq.cell_fluxBC(k)) THEN
                  k1=k
               ENDIF
               IF(kkg.eq.cell_fluxBC(k)) THEN
                  k2=k
               ENDIF               
            ENDDO
            IF(k1.ne.0.and.k2.ne.0) THEN
               IF((zone_fluxBC(k1).eq.ivalve.and.zone_fluxBC(k2).eq.ovalve).or. &
                  (zone_fluxBC(k2).eq.ivalve.and.zone_fluxBC(k1).eq.ovalve)) THEN
!
!                 count number of valves in each MPI domain                      
                  IF(initial_loc) THEN
                     num_valveloc=num_valveloc+1  
                     initial_loc=.false.
                  ENDIF   
!                     
                  IF(num_valveloc.ne.0) THEN
!                    count number of faces of a valve in each local subdomain                                       
                     num_valveface(num_valveloc)=num_valveface(num_valveloc)+1
!
!                    find the face index of a valve
                     nn=num_valveface(num_valveloc)
                     n_face_valve(num_valveloc,nn)=i1                 !non-face index for the valve faces                         
                  ENDIF
!
                  IF(zone_fluxBC(k1).eq.ivalve.and.zone_fluxBC(k2).eq.ovalve) THEN
                     icell_valve(num_valveloc,nn)=ii               !fzone_valve(tt,1) local ii
                     ocell_valve(num_valveloc,nn)=kk               !fzone_valve(tt,2) local kk (kk=1~ncell_fp)
                  ELSEIF(zone_fluxBC(k2).eq.ivalve.and.zone_fluxBC(k1).eq.ovalve) THEN
                     icell_valve(num_valveloc,nn)=kk               !fzone_valve(tt,1) local kk (kk=1~ncell_fp)
                     ocell_valve(num_valveloc,nn)=ii               !fzone_valve(tt,2) local ii
                  ENDIF                      
!                     
               ENDIF
            ENDIF
         ENDDO
      ENDDO  
!
!......Mapping valve faces to a valve
!    
      mapping_valve=0
      DO zz=1,num_valveloc                     !local subdomain valve numbering
         DO i=1,num_valveface(zz)              !local subdomain valve face numbering
            ii=icell_valve(zz,i)
            kk=ocell_valve(zz,i)
            iig=get_global_cell(ii) 
            kkg=get_global_cell(kk)            
!
!             Searching valve number in global domain             
            DO tt=1,num_valve                  !global domain loop                  
               IF(nzone_fluxBC(iig).eq.fzone_valve(tt,1).or.nzone_fluxBC(kkg).eq.fzone_valve(tt,1)) THEN !.fzone_valve(tt,1)=valve
                  mapping_valve(zz)=tt         !local-to-global mapping
               ENDIF
            ENDDO
!
         ENDDO 
      ENDDO    
!      
!.....Find nonk_valve for asymmetric sumation
!   
      nf_number=0
      istart=istart_nf(1,nf_number)  
      DO zz=1,num_valveloc                  !local loop
         tt=mapping_valve(zz)               !mapping local MCP to global MCP number
         DO i=1,num_valveface(zz)
            i1=n_face_valve(zz,i) 
            ir=i1-istart    !i1=istart+i
            ii=left_nf(i1)
            kk=right_non(ir)
!            
            ii_valve=icell_valve(zz,i)
            kk_valve=ocell_valve(zz,i)
!            
            DO j=1,nf_nonk
               k=right_nb_k(j) 
               ii_tmp=left_nf(k)               
               kk_tmp=right_non(k) 
               
               IF((ii_valve.eq.ii_tmp.and.kk_valve.eq.kk_tmp).or.(ii_valve.eq.kk_tmp.and.kk_valve.eq.ii_tmp)) THEN
                  nonk_valve(zz,i)=j  
               ENDIF
            ENDDO              
         ENDDO
      ENDDO      
!      
      DEALLOCATE(cell_fluxBC,zone_fluxBC,nzone_fluxBC)
!      
      RETURN
      END SUBROUTINE  valve_define_fluxBC
