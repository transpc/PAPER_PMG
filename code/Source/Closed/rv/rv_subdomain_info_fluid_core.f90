!
      SUBROUTINE rv_subdomain_info_fluid_core(nelem)
!
!     This routine defines local arrays for each subdomain.  All the 
!     arrays can be automatically defined based on the coloring (manually
!     or METIS). It includes Window METIS (Ver.5.0) with having been 
!     activated in preprocessor option (metis=.true.)
!
      USE Zrv_mpi   , ONLY: celem_fluid_core,jperm_fluid_core, &
                            jjperm_fluid_core,ncell_fluid1_core,ncell_fluid1_core_dsp
      USE Zrv_ncell , ONLY: ncell_fluid_core
      USE Zcore     , ONLY: np,myrank
!
      IMPLICIT NONE      
!
!     input
      INTEGER nelem
!     local variables
      INTEGER :: j,ip,jp,ie
      INTEGER :: cinter
!     local arrays
      INTEGER :: ja(nelem)
!
      IF(np.eq.1) celem_fluid_core(:)=1
!
!     count and collect
!
      ip=myrank+1
      cinter=0
      DO ie=1,nelem
         jp=celem_fluid_core(ie)
         IF(jp.eq.ip)THEN
            cinter=cinter+1
            ja(cinter)=ie
         ENDIF
      ENDDO
!
!     number
!
      ncell_fluid_core=cinter
      ALLOCATE(jperm_fluid_core(ncell_fluid_core))
      DO j=1,cinter
         jperm_fluid_core(j)=ja(j)
      ENDDO
!
!     get ncell_fluid for all procs
!
      ALLOCATE(ncell_fluid1_core(np),ncell_fluid1_core_dsp(np))
      DO ip=1,np
         ncell_fluid1_core(ip)=0
      ENDDO
      DO ie=1,nelem
         ip=celem_fluid_core(ie)
         ncell_fluid1_core(ip)=ncell_fluid1_core(ip)+1
      ENDDO
!
!     get size and displacement to call ALLGATHERV
!
         ip=1
         ncell_fluid1_core_dsp(ip)=0
      DO ip=2,np
         ncell_fluid1_core_dsp(ip)=ncell_fluid1_core_dsp(ip-1)+ncell_fluid1_core(ip-1)
      ENDDO
!
!     build jjperm out of local jperm
!
      ALLOCATE(jjperm_fluid_core(nelem))
      IF(np.gt.1) THEN
         CALL allgather_vec_i(jperm_fluid_core,ncell_fluid_core,jjperm_fluid_core,nelem,ncell_fluid1_core,ncell_fluid1_core_dsp)
      ELSE
         DO j=1,nelem
            jjperm_fluid_core(j)=jperm_fluid_core(j)
         ENDDO
      ENDIF
!
      END SUBROUTINE rv_subdomain_info_fluid_core
!
