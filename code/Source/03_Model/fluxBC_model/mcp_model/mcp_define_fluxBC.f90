!
      SUBROUTINE mcp_define_fluxBC
!
!     DATA STRUCTURE
!
!     num_mcp=total number of MCPs in global domain
!     num_mcploc=number of MCPs in local subdomain
!     num_mcpface(1:num_mcploc)=number of faces of an MCP in local subdomain
!     n_face_mcp(1:num_mcploc,1:num_mcpface)=face index of an MCP in local subdomain
!     fzone_mcp(1:num_mcp,1)=fluid zone index of MCP in global domain (rv_parameters.in)
!     fzone_mcp(1:num_mcp,2)=fluid zone index of MCP in global domain (rv_parameters.in)
!     icell_mcp(1:num_mcploc,1:num_mcpface)=local domain cell index of an MCP
!     ocell_mcp(1:num_mcploc,1:num_mcpface)=local domain cell index of an MCP
!     mcp_area(1:num_mcploc)=sum of area of each MCP face in local subdomain
!     mcp_on(1:num_mcp): wheather mcp is in operation(1) or not (0) (global domain)
!     dir_face_mcp(1:num_mcploc,1:num_mcpface)=local domain cell index of an MCP
!  
      USE Zinterface
      USE Zzone        , ONLY: nzone,ncell_fluid_all
      USE Zcore        , ONLY: np,myrank
      USE Znum_cell    , ONLY: istart_nf,right_nb_k
      USE Zvec_index   , ONLY: left_nf,right_non 
      USE Zvec_param   , ONLY: nf_nonk
      USE Zparam       , ONLY: pi
      USE Zmcp         , ONLY: nonk_mcp,num_mcp,num_mcpface,num_mcploc,icell_mcp,dir_face_mcp,mapping_mcp, &
                               ocell_mcp,n_face_mcp,fluxd_mcpface,fluxg_mcpface,fluxl_mcpface,  &
                               fzone_mcp
!
      IMPLICIT NONE
!
      INTEGER :: get_global_cell
!      
      INTEGER :: nf_number,istart,isize,i1
      INTEGER :: i,ii,kk,tt,nn,zz,k,ir,j,ii_tmp,kk_tmp
      INTEGER :: imcp,omcp,ii_mcp,kk_mcp
      INTEGER :: err,unit_valve,num_zone,iig,kkg,k1,k2
      INTEGER,DIMENSION(:),ALLOCATABLE :: cell_fluxBC,zone_fluxBC,nzone_fluxBC,nn_array
      LOGICAL,SAVE :: initial_loc
      LOGICAL,SAVE :: initial_nzone=.true.
!
!.....Open FluxBCzone file
!      
      unit_valve=424 
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
!......Define global nzone
!      
      IF(initial_nzone) THEN
         PRINT*,'>>>MCP face is finding...'
         IF(np.gt.1) CALL communicate_1d_int(nzone)         
         initial_nzone=.false.
      ENDIF   
!
      num_mcploc=0
      num_mcpface=0
      nf_number=0
      istart=istart_nf(1,nf_number)
      isize =istart_nf(2,nf_number)
      DO tt=1,num_mcp             !num_mcp=total number of mcps of a reactor        
         imcp=fzone_mcp(tt,1)     !inflow fluid zone number (MCP)
         omcp=fzone_mcp(tt,2)     !outflow fluid zone number (Distributor)
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
               IF((zone_fluxBC(k1).eq.imcp.and.zone_fluxBC(k2).eq.omcp).or. &
                  (zone_fluxBC(k2).eq.imcp.and.zone_fluxBC(k1).eq.omcp)) THEN            
!                   
                  IF(initial_loc) THEN
                     num_mcploc=num_mcploc+1  
                     initial_loc=.false.
                  ENDIF  
                  IF(num_mcploc.ne.0) num_mcpface(num_mcploc)=num_mcpface(num_mcploc)+1               
               ENDIF
            ENDIF          
         ENDDO     
      ENDDO         
!
!....Initialize the array wrt the max value of mcp faces
!      
      ALLOCATE(nn_array(np))
      nn_array=0
      nn_array(myrank+1)=MAXVAL(num_mcpface(1:num_mcp)) 
      CALL allreducei_i(nn_array,np)       
      nn=MAXVAL(nn_array)      
!      nn=MAXVAL(num_mcpface(1:num_mcp))
      IF(.not.ALLOCATED(n_face_mcp)) ALLOCATE(n_face_mcp(num_mcp,nn))     !local domain 
      IF(.not.ALLOCATED(icell_mcp)) ALLOCATE(icell_mcp(num_mcp,nn))       !local domain
      IF(.not.ALLOCATED(ocell_mcp)) ALLOCATE(ocell_mcp(num_mcp,nn))       !local domain 
      IF(.not.ALLOCATED(fluxl_mcpface)) ALLOCATE(fluxl_mcpface(num_mcp))  !global domain, initialization is in mcp_flow.f90
      IF(.not.ALLOCATED(fluxg_mcpface)) ALLOCATE(fluxg_mcpface(num_mcp))  !global domain, initialization is in mcp_flow.f90
      IF(.not.ALLOCATED(fluxd_mcpface)) ALLOCATE(fluxd_mcpface(num_mcp))  !global domain, initialization is in mcp_flow.f90
      IF(.not.ALLOCATED(mapping_mcp)) ALLOCATE(mapping_mcp(num_mcp))      !local-to-global
      IF(.not.ALLOCATED(dir_face_mcp)) ALLOCATE(dir_face_mcp(num_mcp,nn)) !local domain
      IF(.not.ALLOCATED(nonk_mcp)) ALLOCATE(nonk_mcp(num_mcp,nn))     !local domain 
!      
      num_mcploc=0
      num_mcpface=0
      n_face_mcp=0
      icell_mcp=0
      ocell_mcp=0
      nonk_mcp=0
!      
!.....Finding mcp faces
!
      nn=0
      istart=istart_nf(1,nf_number)
      isize =istart_nf(2,nf_number)          
      DO tt=1,num_mcp
         imcp=fzone_mcp(tt,1)     !inflow fluid zone number of an mcp
         omcp=fzone_mcp(tt,2)     !outflow fluid zone number of an mcp
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
               IF((zone_fluxBC(k1).eq.imcp.and.zone_fluxBC(k2).eq.omcp).or. &
                  (zone_fluxBC(k2).eq.imcp.and.zone_fluxBC(k1).eq.omcp)) THEN
!
!                 count number of valves in each MPI domain                      
                  IF(initial_loc) THEN
                     num_mcploc=num_mcploc+1   
                     initial_loc=.false.
                  ENDIF   
!                     
                  IF(num_mcploc.ne.0) THEN
!                    count number of faces of an mcp in each local subdomain                                       
                     num_mcpface(num_mcploc)=num_mcpface(num_mcploc)+1
!
!                    find the face index of an mcp
                     nn=num_mcpface(num_mcploc)
                     n_face_mcp(num_mcploc,nn)=i1                 !non-face index for the mcp faces                         
                  ENDIF
!
                  IF(zone_fluxBC(k1).eq.imcp.and.zone_fluxBC(k2).eq.omcp) THEN
                     icell_mcp(num_mcploc,nn)=ii               !fzone_mcp(tt,1) local ii
                     ocell_mcp(num_mcploc,nn)=kk               !fzone_mcp(tt,2) local kk (kk=1~ncell_fp)
                  ELSEIF(zone_fluxBC(k2).eq.imcp.and.zone_fluxBC(k1).eq.omcp) THEN
                     icell_mcp(num_mcploc,nn)=kk               !fzone_mcp(tt,1) local kk (kk=1~ncell_fp)
                     ocell_mcp(num_mcploc,nn)=ii               !fzone_mcp(tt,2) local ii
                  ENDIF      
!                  
                  IF(zone_fluxBC(k1).lt.zone_fluxBC(k2)) THEN
                     dir_face_mcp(num_mcploc,nn)=1.d0    !flux direction is in increasing order of fluid zone
                  ELSE
                     dir_face_mcp(num_mcploc,nn)=-1.d0   !correction of flux direction to make MCP flux in in increasing order of fluid zone
                  ENDIF                   
!                     
               ENDIF
            ENDIF
         ENDDO
      ENDDO  
!
!.....mapping local MCPs to global MCP number
!      
      mapping_mcp=0
      DO zz=1,num_mcploc                     !local subdomain valve numbering
         DO i=1,num_mcpface(zz)              !local subdomain valve face numbering
            ii=icell_mcp(zz,i)
            kk=ocell_mcp(zz,i)
            iig=get_global_cell(ii) 
            kkg=get_global_cell(kk)            
!
!             Searching valve number in global domain             
            DO tt=1,num_mcp                  !global domain loop                  
               IF(nzone_fluxBC(iig).eq.fzone_mcp(tt,1).or.nzone_fluxBC(kkg).eq.fzone_mcp(tt,1)) THEN !.fzone_valve(tt,1)=valve
                  mapping_mcp(zz)=tt         !local-to-global mapping
               ENDIF
            ENDDO
!
         ENDDO 
      ENDDO        
!      
!.....Find nonk_mcp for asymmetric sumation
!         
      nf_number=0
      istart=istart_nf(1,nf_number)  
      DO zz=1,num_mcploc                  !local loop
         tt=mapping_mcp(zz)               !mapping local MCP to global MCP number
         DO i=1,num_mcpface(zz)
            i1=n_face_mcp(zz,i) 
            ir=i1-istart    !i1=istart+i
            ii=left_nf(i1)
            kk=right_non(ir)
!            
            ii_mcp=icell_mcp(zz,i)
            kk_mcp=ocell_mcp(zz,i)
!            
            DO j=1,nf_nonk
               k=right_nb_k(j) 
               ii_tmp=left_nf(k)               
               kk_tmp=right_non(k) 
               
               IF((ii_mcp.eq.ii_tmp.and.kk_mcp.eq.kk_tmp).or.(ii_mcp.eq.kk_tmp.and.kk_mcp.eq.ii_tmp)) THEN
                  nonk_mcp(zz,i)=j  
               ENDIF
   
            ENDDO              
         ENDDO
      ENDDO
      
      
      DEALLOCATE(cell_fluxBC,zone_fluxBC,nzone_fluxBC)      
!      
      END SUBROUTINE  mcp_define_fluxBC
