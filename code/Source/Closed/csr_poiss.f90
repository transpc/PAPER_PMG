!
      SUBROUTINE csr_build_a_c(poiss_diag,poiss_non_i,poiss_non_k)
!
      USE Zzone         , ONLY: ncell_cond
      USE Zvec_param    , ONLY: nfc_nonk,nfc_non
      USE Znum_cell     , ONLY: iac_nb,iptrc_nb_k
      USE Zmpi          , ONLY: au_c,ju_a_c,ia_a_c
!     USE Zvec_index_solid    , ONLY: id_nf,right_solid_nf
!
      IMPLICIT NONE
!.....Input
      REAL(8),DIMENSION(ncell_cond) :: poiss_diag
      REAL(8),DIMENSION(nfc_non) :: poiss_non_i
      REAL(8),DIMENSION(nfc_nonk) :: poiss_non_k
!.....Local variables
      INTEGER i,i1,j
!     INTEGER k1
      INTEGER ip1,ip2
!
      ip1=1
      DO i=1,ncell_cond
!........Lower part, index to ia_a set to zero if no lower part
         i1=iptrc_nb_k(i)
         IF(i1.gt.0) then
            ip2=ia_a_c(i)
            DO j=iac_nb(i1),iac_nb(i1+1)-1
               au_c(ip2)=poiss_non_k(j)
!              k1=left_nf(k)
!              if(npb(k1).gt.0) au(ip2)=0.0d0
               ip2=ip2+1
            ENDDO
         ENDIF
!........Diagonal
         j=ju_a_c(i)
         au_c(j)=poiss_diag(i)
         poiss_diag(i)=1.d0/poiss_diag(i)
!........Upper part
         DO j=ju_a_c(i)+1,ia_a_c(i+1)-1
            au_c(j)=poiss_non_i(ip1)
!           k=right_solid_nf(ip1)
!           if(npb(k).gt.0) au(j)=0.0d0
            ip1=ip1+1
         ENDDO
      ENDDO
!
      END SUBROUTINE csr_build_a_c
!
      SUBROUTINE csr_build_a(poiss_diag,poiss_non_i,poiss_non_k)
!
      USE Zzone         , ONLY: ncell_fluid
      USE Zvec_param    , ONLY: nf_nonk,nf_non
      USE Znum_cell     , ONLY: ia_nb,iptr_nb_k
      USE Zmpi          , ONLY: au,ju_a,ia_a
!     USE Zvec_index    , ONLY: left_nf,right_non
!     USE Zbc_index     , ONLY: npb
!
      IMPLICIT NONE
!     input
      REAL(8),DIMENSION(ncell_fluid) :: poiss_diag
      REAL(8),DIMENSION(nf_non) :: poiss_non_i
      REAL(8),DIMENSION(nf_nonk) :: poiss_non_k
!     local variables
      INTEGER i,i1,j
!     INTEGER k1
      INTEGER ip1,ip2
!
      ip1=1
      DO i=1,ncell_fluid
!........Lower part, index to ia_a set to zero if no lower part
         i1=iptr_nb_k(i)
         IF(i1.gt.0) then
            ip2=ia_a(i)
            DO j=ia_nb(i1),ia_nb(i1+1)-1
               au(ip2)=poiss_non_k(j)
!              k1=left_nf(k)
!              if(npb(k1).gt.0) au(ip2)=0.0d0
               ip2=ip2+1
            ENDDO
         ENDIF
!........Diagonal
         j=ju_a(i)
         au(j)=poiss_diag(i)
         poiss_diag(i)=1.d0/poiss_diag(i)
!........Upper part
         DO j=ju_a(i)+1,ia_a(i+1)-1
            au(j)=poiss_non_i(ip1)
!           k=right_non(ip1)
!           if(npb(k).gt.0) au(j)=0.0d0
            ip1=ip1+1
         ENDDO
      ENDDO
!
      END SUBROUTINE csr_build_a
