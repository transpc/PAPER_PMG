!
      SUBROUTINE udfn_s_gapij
!
!     Modifies the momentum source terms at the free surface cells
!
      USE Zporous    , ONLY: sgap,s_gapij_non_i,s_gapij_non_k
!
      USE Znum_cell  , ONLY: istart_nf
      USE Zvec_index , ONLY: left_nf,right_non
      USE Zvec_geo   , ONLY: xn_nf
!
      IMPLICIT NONE
!            
      INTEGER :: i
      INTEGER :: ii,kk
      INTEGER :: nf_number,istart,len,i1
!      
!====> non only
      nf_number=0
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)

      DO i=1,len
         i1=istart+i
         ii=left_nf(i1)
         kk=right_non(i)
         IF(abs(xn_nf(i1,1)).eq.1.0d0)THEN
!====> ii
            s_gapij_non_i(i)=sgap(ii,1)
!====> kk
            s_gapij_non_k(i)=sgap(kk,1)
         ELSEIF(abs(xn_nf(i1,2)).eq.1.0d0)THEN
!====> ii
            s_gapij_non_i(i)=sgap(ii,2)
!====> kk
            s_gapij_non_k(i)=sgap(kk,2)
         ENDIF
      ENDDO
!
      RETURN
      END SUBROUTINE udfn_s_gapij
