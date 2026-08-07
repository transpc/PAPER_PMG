!
      SUBROUTINE subdomain_info_ser(nelem) 
!
!     This routine defines local arrays for each subdomain.  All the 
!     arrays can be automatically defined based on the coloring (manually
!     or METIS). It includes Window METIS (Ver.5.0) with having been 
!     activated in preprocessor option (metis=.true.)
!
      USE Zinterface
      USE Zmpi        , ONLY: celem,ncell_fp,                               &
                              maxmt_fp,maxmt_fluid,maxmt_cell,maxmt_nfluid, &
                              jperm,jjperm,                                 &
                              i3perm
      USE Zzone       , ONLY: ncell_fluid
      USE Zparam      , ONLY: ns,nin_max,nb_max
      USE Zbc_index   , ONLY: nbcon
      USE Znum_cell   , ONLY: i_neigh_tmp,j_neigh_tmp,j_nbcon_tmp,index_sort_tmp1, &
                              i_neigh,neigh,index_sort,indexr_sort
      USE Znormal     , ONLY: nji,i_neigh_nbcon0,maxmt_fluid_nbcon0
! PMG 
      USE MD_parameter, ONLY: nf_max
      USE MD_geometry, ONLY: nelem_mg,num_neigh_mg,neigh_mg
!
      IMPLICIT NONE      
!
!.....Input
      INTEGER :: nelem
!.....Local variables
      INTEGER :: i,j,k,jj,ne
      INTEGER :: j0,j1,j2,k0,l
!.....Local arrays
      INTEGER :: neighl(ns)
      INTEGER :: temp_sort(ns)
!
      celem=1
      ncell_fluid=nelem
      ncell_fp=ncell_fluid
      ALLOCATE(jperm(ncell_fluid))
      ALLOCATE(jjperm(nelem))
      DO jj=1,nelem
         jperm(jj)=jj
         jjperm(jj)=jj
      ENDDO
!
      ALLOCATE(i_neigh(ncell_fluid+1))
      i_neigh(1)=1
      DO i=1,ncell_fluid
         i_neigh(i+1)=i_neigh_tmp(i+1)
      ENDDO
      maxmt_fluid=i_neigh(ncell_fluid+1)-1
      ALLOCATE(neigh(maxmt_fluid))
      ALLOCATE(nbcon(maxmt_fluid))
      ALLOCATE(index_sort(maxmt_fluid),indexr_sort(maxmt_fluid))
      DO i=1,ncell_fluid
         j1=i_neigh(i)-1
         l=i_neigh(i+1)-i_neigh(i)
!
!........Sort in the vector order
!
         DO j=1,l
            index_sort(j+j1)=j
            IF(j_nbcon_tmp(j1+j).eq.0)THEN
               k=j_neigh_tmp(j1+j)
               neighl(j)=k
            ELSE
               IF    (j_nbcon_tmp(j1+j).ge.201)THEN
                  neighl(j)=ncell_fluid+1
               ELSEIF(j_nbcon_tmp(j1+j).gt.0.and.j_nbcon_tmp(j1+j).le.nin_max)THEN
                  neighl(j)=ncell_fluid+2
               ELSEIF(j_nbcon_tmp(j1+j).gt.nin_max.and.j_nbcon_tmp(j1+j).le.nb_max)THEN
                  neighl(j)=ncell_fluid+3
               ELSEIF(j_nbcon_tmp(j1+j).eq.-1)THEN
                  neighl(j)=ncell_fluid+4
               ELSEIF(j_nbcon_tmp(j1+j).eq.-2)THEN
                  neighl(j)=ncell_fluid+5
               ELSEIF(j_nbcon_tmp(j1+j).eq.-3.or.j_nbcon_tmp(j1+j).eq.-4)THEN
                  neighl(j)=ncell_fluid+6
               ELSEIF(j_nbcon_tmp(j1+j).eq.-5.or.j_nbcon_tmp(j1+j).eq.-6)THEN
                  neighl(j)=ncell_fluid+7
               ELSEIF(j_nbcon_tmp(j1+j).eq.101)THEN
                  neighl(j)=ncell_fluid+8
               ENDIF
            ENDIF
         ENDDO
         CALL sortx_i(neighl,index_sort(j1+1),l)
!
!........Restore neigh to 0 for non compute
!
         DO j=l,1,-1
            if(neighl(j).le.ncell_fluid) exit
            neighl(j)=0
         ENDDO
!........Get sorted nbcon
         DO j=1,l
            k=index_sort(j+j1)
            indexr_sort(k+j1)=j
            temp_sort(j)=j_nbcon_tmp(j1+k)
         ENDDO
         DO j=1,l
            neigh(j+j1)=neighl(j)
            nbcon(j+j1)=temp_sort(j)
         ENDDO
      ENDDO
!
      ALLOCATE(i_neigh_nbcon0(ncell_fluid+1))
      i_neigh_nbcon0(1)=1
      j0=1
      DO i=1,ncell_fluid
         DO j=i_neigh(i),i_neigh(i+1)-1
            IF(nbcon(j).ne.0) EXIT
            j0=j0+1
         ENDDO
         i_neigh_nbcon0(i+1)=j0
      ENDDO
      maxmt_fluid_nbcon0=i_neigh_nbcon0(ncell_fluid+1)-1
!
      ALLOCATE(nji(maxmt_fluid_nbcon0))
      DO i=1,ncell_fluid
         j2=i_neigh_nbcon0(i)
         DO j=i_neigh(i),i_neigh(i+1)-1
            IF(nbcon(j).ne.0) EXIT
            k=neigh(j)
            k0=i_neigh(k)-1
            DO j1=i_neigh(k),i_neigh(k+1)-1
               IF(neigh(j1).eq.i) then
                   nji(j2)=j1-k0
                   j2=j2+1
                   EXIT
                ENDIF
            ENDDO
         ENDDO
      ENDDO
!
!.....Get index_sort global
!
      ALLOCATE(index_sort_tmp1(maxmt_cell))
      DO i=1,maxmt_nfluid
         index_sort_tmp1(i)=index_sort(i)
      ENDDO
!
!.....Get sorted global j_neigh_tmp to get nji
!.....Get sorted global j_nbcon_tmp to be used in solid
!
      DO ne=1,nelem
         j0=i_neigh_tmp(ne)-1
         DO j=i_neigh_tmp(ne),i_neigh_tmp(ne+1)-1
            k=index_sort_tmp1(j)
            temp_sort(j-j0)=j_neigh_tmp(j0+k)
         ENDDO
         DO j=i_neigh_tmp(ne),i_neigh_tmp(ne+1)-1
            j_neigh_tmp(j)=temp_sort(j-j0)
         ENDDO
         DO j=i_neigh_tmp(ne),i_neigh_tmp(ne+1)-1
            k=index_sort_tmp1(j)
            temp_sort(j-j0)=j_nbcon_tmp(j0+k)
         ENDDO
         DO j=i_neigh_tmp(ne),i_neigh_tmp(ne+1)-1
            j_nbcon_tmp(j)=temp_sort(j-j0)
         ENDDO
      ENDDO
!
      maxmt_fp=i_neigh(ncell_fluid+1)-1
!
!     Global cell to local cell for CUPID_RV
!
      ALLOCATE(i3perm(nelem))
      i3perm(:)=0
      DO i=1,ncell_fluid
         i3perm(jperm(i))=i
      ENDDO     
      
! PMG copy to PMG solver
!      ndom = np

      nelem_mg = nelem
      nf_max = ns
      allocate(num_neigh_mg(nelem_mg),neigh_mg(nf_max,nelem_mg))

      DO i=1,nelem
         k = 0
         DO j=i_neigh(i),i_neigh(i+1)-1
             j0 = neigh(j)
            IF(j0.EQ.0) EXIT
            k = k+1
            neigh_mg(k,i) = j0
         ENDDO 
         num_neigh_mg(i) = k
      ENDDO
! - - - - - - - - - - - - - - !
!
!     Deallocate the local variables for making subdomain array
!
!     DEALLOCATE(irecv,jsend)
!
      END SUBROUTINE subdomain_info_ser
