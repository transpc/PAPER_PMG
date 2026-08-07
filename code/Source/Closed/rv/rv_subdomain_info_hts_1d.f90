!
      SUBROUTINE rv_subdomain_info_hts_1d(nelem)
!
!     This routine defines local arrays for each subdomain.  All the 
!     arrays can be automatically defined based on the coloring (manually
!     or METIS). It includes Window METIS (Ver.5.0) with having been 
!     activated in preprocessor option (metis=.true.)
!
      USE Zrv_mpi    , ONLY: celem_hts_1d,jperm_hts_1d
      USE Zrv_hts_1d , ONLY: ncell_hts_1d
      USE Zcore      , ONLY: np,myrank
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
      IF(np.eq.1) celem_hts_1d(:)=1
!
!     count and collect
!
      ip=myrank+1
      cinter=0
      DO ie=1,nelem
         jp=celem_hts_1d(ie)
         IF(jp.eq.ip)THEN
            cinter=cinter+1
            ja(cinter)=ie
         ENDIF
      ENDDO
!
!     number
!
      ncell_hts_1d=cinter
      ALLOCATE(jperm_hts_1d(ncell_hts_1d))
      DO j=1,cinter
         jperm_hts_1d(j)=ja(j)
      ENDDO
!
      END SUBROUTINE rv_subdomain_info_hts_1d
!
