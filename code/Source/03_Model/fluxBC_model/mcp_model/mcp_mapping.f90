!
      SUBROUTINE mcp_mapping
!
      USE Zcore        , ONLY: np
      USE Zvec_geo     , ONLY: saa_nf
      USE Zcoord3      , ONLY: volp
      USE Zmcp         
      USE Zvec_index   , ONLY: left_nf,right_non
      USE Znum_cell    , ONLY: istart_nf
      USE Zmpi         , ONLY: celem
!
      IMPLICIT NONE
!
!......External function
      INTEGER :: get_global_cell      
      
      INTEGER :: tt,i,ii,kk,zz,i1,ii1,kk1,ir,istart,nf_number
!      
!.....Calculate the mcp area and volume (local subdomain to global domain)
!      
      IF(.not.ALLOCATED(mcp_area)) ALLOCATE(mcp_area(num_mcp)) !num_mcp is always bigger than num_mcploc
      IF(.not.ALLOCATED(mcp_vol)) ALLOCATE(mcp_vol(num_mcp))   !num_mcp is always bigger than num_mcploc     
      IF(.not.ALLOCATED(mcp_area_global)) ALLOCATE(mcp_area_global(num_mcp)) !num_mcp is always bigger than num_mcploc
      IF(.not.ALLOCATED(mcp_vol_global)) ALLOCATE(mcp_vol_global(num_mcp))   !num_mcp is always bigger than num_mcploc    
      IF(.not.ALLOCATED(num_mcpface_global)) ALLOCATE(num_mcpface_global(num_mcp))   !num_mcp is always bigger than num_mcploc     
      IF(.not.ALLOCATED(num_mcpface_tmp)) ALLOCATE(num_mcpface_tmp(num_mcp))   !num_mcp is always bigger than num_mcploc     
      mcp_area(:)=0.d0
      mcp_vol(:)=0.d0
      num_mcpface_tmp(:)=0
      mcp_area_global(:)=0.d0
      mcp_vol_global(:)=0.d0  
      num_mcpface_global(:)=0
      
      nf_number=0
      istart=istart_nf(1,nf_number)   
!      
      DO zz=1,num_mcploc
         tt=mapping_mcp(zz)
         DO i=1,num_mcpface(zz)
            i1=n_face_mcp(zz,i)
            ir=i1-istart
 !           
            ii=left_nf(i1)
            kk=right_non(ir)
            ii1=get_global_cell(ii)
            kk1=get_global_cell(kk)
!
            IF(celem(ii1).ne.celem(kk1)) THEN  !MPI Processors for MCP and discharger are different.
               mcp_area(tt)=mcp_area(tt)+saa_nf(i1)*0.5d0
               mcp_vol(tt)=mcp_vol(tt)+volp(ii)*0.5d0
               num_mcpface_tmp(tt)=num_mcpface_tmp(tt)+0.5d0
            ELSE
               mcp_area(tt)=mcp_area(tt)+saa_nf(i1)
               mcp_vol(tt)=mcp_vol(tt)+volp(ii)                
               num_mcpface_tmp(tt)=num_mcpface_tmp(tt)+1.d0
            ENDIF
         ENDDO
      ENDDO  
!   
      IF(np.gt.1) THEN
         CALL allreduce_r(mcp_area,mcp_area_global,num_mcp)
         CALL allreduce_r(mcp_vol,mcp_vol_global,num_mcp)
         CALL allreduce_r(num_mcpface_tmp,num_mcpface_global,num_mcp)
      ELSE
         mcp_area_global=mcp_area 
         mcp_vol_global=mcp_vol
         num_mcpface_global=num_mcpface
      ENDIF  
!
      RETURN
      END SUBROUTINE mcp_mapping