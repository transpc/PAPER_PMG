!
      SUBROUTINE calc_solid
! 
!     This routine calculates solid temperature in porous medium. 
!
      USE SOLID_DATA    , ONLY: solid
      USE Zmpi          , ONLY: ncell_ps,maxmt_c,maxmt_pad_c,                          &
                                max_neigh_c,                                           &
                                au_c,ju_a_c,ja_a_c,ia_a_c,iend_c,                      &
                                ap_c,iap_c,jap_c,jaa_c,jaar_c,iaa_c,ngroup_c,nbgroup_c
      USE Zvec_param    , ONLY: nfc_nonk,nfc_non
      USE Zzone         , ONLY: ncell_cond,ncell_cond_pad,ncell_cond_padv,nmaterial_c
      USE Zbc_index     , ONLY: nbcon_c
      USE Zbicg         , ONLY: pbcgind
      USE Znum_cell     , ONLY: neigh_c,i_neigh_c
      USE Zuserdefined  , ONLY: udfl_mat_prop
      USE Zcore         , ONLY: myrank
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i,j,j0
      INTEGER :: iOKr,iOKk
      LOGICAL, SAVE:: INITIAL=.true.
      REAL(8) :: CpVol,Condu
      REAL(8) :: epsFactor= 1.d0
!.....Local array
      REAL(8),DIMENSION(ncell_cond) :: poiss_diag_solid,src_solid
      REAL(8),DIMENSION(ncell_ps) :: v1p_c
!.....Local vector array
      REAL(8),DIMENSION(nfc_non) :: poiss_solid_non_i
      REAL(8),DIMENSION(nfc_nonk) :: poiss_solid_non_k

!      
!.....Revise solid properties such as density*specific heat and conductivty
!
      IF(udfl_mat_prop)THEN
         CALL udfn_mat_prop
      ELSE 
          DO i=1,ncell_cond
             IF(IABS(nmaterial_c(i)) .lt. 50) THEN
                CALL mat_prop(IABS(nmaterial_c(i)),solid%tsol(i),CpVol,Condu,iOKr,iOKk)
                solid%rhocps(i)=CpVol 
                solid%conds(i) =Condu 
             ENDIF
             solid%matnum(i)=IABS(nmaterial_c(i))
          ENDDO      
      ENDIF
!
!.....Calculate solid temperature
!.....Set matrix for solid temperature equation
!
      CALL scalar_matrix_solid(poiss_diag_solid,poiss_solid_non_i,poiss_solid_non_k,src_solid)
!
!........Set CSR format for neigh(:,:). neigh CSR is defined only ONCE.
!
      IF(INITIAL)THEN
         INITIAL=.FALSE.
!
!........Get the maxmt to allocate ja_a
!
         Allocate(ia_a_c(ncell_cond+1),ju_a_c(ncell_cond),iend_c(ncell_cond))
         ia_a_c(1)=1
         maxmt_c=0
         max_neigh_c=0
         DO i=1,ncell_cond
            j0=i_neigh_c(i)-1
!
!...........Get neighbors count for compute elements only
!
            DO j=i_neigh_c(i),i_neigh_c(i+1)-1
               IF(nbcon_c(j).ne.0) exit
               maxmt_c=maxmt_c+1
            ENDDO
!
!...........Add 1 for new diagonal entry
!
            maxmt_c=maxmt_c+1
            ia_a_c(i+1)=maxmt_c+1
            max_neigh_c=max(max_neigh_c,ia_a_c(i+1)-ia_a_c(i))
         ENDDO
!
         Allocate(ja_a_c(maxmt_c))
         CALL csr_neigh_c(maxmt_c,ncell_cond,ia_a_c,ja_a_c,ju_a_c,iend_c,i_neigh_c,neigh_c,nbcon_c)
!
!........Process ilup 
!
         CALL reorder_ilup_c
!
!....... Rearrange A in blocks of equal number of neighbors to vectorize
!        write(*,*) 'ngroup_c',ncell_cond,ngroup_c,max_neigh_c
         CALL gener_vect_size(ncell_cond,max_neigh_c,ia_a_c,ngroup_c)
         ALLOCATE(au_c(maxmt_c))
         ALLOCATE(iap_c(2,ngroup_c+1))
         ALLOCATE(iaa_c(2,ngroup_c+1),jaa_c(ncell_cond),jaar_c(ncell_cond))
         ALLOCATE(nbgroup_c(3,ngroup_c))
         CALL gener_vect_u(ncell_cond,maxmt_pad_c,ncell_cond_pad,ncell_cond_padv,max_neigh_c, &
                           ia_a_c,iaa_c,iap_c,                                                &
                           jaa_c,jaar_c,ngroup_c,nbgroup_c)
         ALLOCATE(jap_c(maxmt_pad_c))
         ALLOCATE(ap_c(maxmt_pad_c))
         CALL copy_ja_vector(ncell_cond,maxmt_c,maxmt_pad_c, &
                             ja_a_c,ia_a_c,jap_c,iap_c,        &
                             jaa_c,iaa_c,ngroup_c,nbgroup_c)
      ENDIF
!
!.....Build directly solverCSR  array here
!
      CALL csr_build_a_c(poiss_diag_solid,poiss_solid_non_i,poiss_solid_non_k)
!
!.....Solve
!
      CALL cupid_solvers_c(epsFactor,poiss_diag_solid,src_solid,v1p_c)
!
      IF(pbcgind.gt.0)THEN
         pbcgind=0
         IF(myrank.eq.0)WRITE(*,*)'Iteration fails in solid!!!'
      ENDIF
! 
!.....Update solid temperature
!
      DO i=1,ncell_ps
         solid%tsol(i)=v1p_c(i)
      ENDDO
!
      END SUBROUTINE calc_solid
