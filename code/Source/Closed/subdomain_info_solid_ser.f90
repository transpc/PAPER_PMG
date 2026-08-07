!
      SUBROUTINE subdomain_info_solid_ser(nelem)
!
!     This routine defines local arrays for each subdomain.  All the 
!     arrays can be automatically defined based on the coloring (manually
!     or METIS). It includes Window METIS (Ver.5.0) with having been 
!     activated in preprocessor option (metis=.true.)
!
      USE Zinterface
      USE Zmpi      , ONLY: celem_c,ncell_ps,     &
                            maxmt_ncond,maxmt_ps, &
                            jperm_c,jjperm_c,     &
                            maxmt_nncond
      USE Zzone     , ONLY: ncell_cond
      USE Zparam    , ONLY: ns
      USE Zbc_index , ONLY: nbcon_c
      USE Znum_cell , ONLY: neigh_c,i_neigh_c_tmp,j_neigh_c_tmp,j_nbcon_c_tmp, &
                            i_neigh_c,index_sort_c
      USE Znormal   , ONLY: nji_c
!
      IMPLICIT NONE      
!
!.....Input
      INTEGER :: nelem
!     local variables
      INTEGER :: i,j,k,ne
      INTEGER :: j0,j1,j2,j3
!.....Local arrays
      INTEGER temp_sort(ns)
      INTEGER :: index_sort0(maxmt_nncond)
!
      celem_c=1
      ncell_cond=nelem
      ncell_ps=ncell_cond
!
      ALLOCATE(jperm_c(nelem))
      ALLOCATE(jjperm_c(nelem))
      DO j=1,nelem
         jperm_c(j)=j
         jjperm_c(j)=j
      ENDDO      
!
!
!.....Define local neigh & num_neigh for fluid region
!.....Define local nji
!
      ALLOCATE(i_neigh_c(ncell_ps+1))
      i_neigh_c(1)=1
      DO i=1,ncell_cond
         ne=i
         i_neigh_c(i+1)=i_neigh_c(i)+(i_neigh_c_tmp(ne+1)-i_neigh_c_tmp(ne))
      ENDDO
      maxmt_ncond=i_neigh_c(ncell_cond+1)-1
      ALLOCATE(index_sort_c(maxmt_ncond))
      ALLOCATE(neigh_c(maxmt_ncond))
      ALLOCATE(nbcon_c(maxmt_ncond))
!
      DO i=1,ncell_cond
         ne=i
         j0=i_neigh_c_tmp(ne)-1
         j1=i_neigh_c(i)-1
!
!........Sort in the vector order
!
         DO j=i_neigh_c(i),i_neigh_c(i+1)-1
            index_sort_c(j)=j-j1
            IF(j_nbcon_c_tmp(j-j1+j0).eq.0)THEN
               k=j_neigh_c_tmp(j-j1+j0)
               neigh_c(j)=k
            ELSE
               IF    (j_nbcon_c_tmp(j-j1+j0).eq.-2)THEN
                  neigh_c(j)=ncell_ps+1
               ELSEIF(j_nbcon_c_tmp(j-j1+j0).eq.-3.or.j_nbcon_c_tmp(j-j1+j0).eq.-4)THEN
                  neigh_c(j)=ncell_ps+2
               ELSEIF(j_nbcon_c_tmp(j-j1+j0).eq.-5.or.j_nbcon_c_tmp(j-j1+j0).eq.-6)THEN
                  neigh_c(j)=ncell_ps+3
               ELSE
                  neigh_c(j)=ncell_ps+4
               ENDIF
            ENDIF
         ENDDO
         CALL sortx_i(neigh_c(j1+1),index_sort_c(j1+1),i_neigh_c(i+1)-i_neigh_c(i))
!
!........Restore neigh to 0 for non compute
!
         DO j=i_neigh_c(i+1)-1,i_neigh_c(i),-1
            if(neigh_c(j).le.ncell_ps) exit
            neigh_c(j)=0
         ENDDO
!........Get sorted nbcon
         DO j=i_neigh_c(i),i_neigh_c(i+1)-1
            k=index_sort_c(j)
            temp_sort(j-j1)=j_nbcon_c_tmp(j0+k)
         ENDDO
         DO j=i_neigh_c(i),i_neigh_c(i+1)-1
            nbcon_c(j)=temp_sort(j-j1)
         ENDDO
      ENDDO
!
!.....Get index_sort global
!
      DO i=1,maxmt_nncond
         index_sort0(i)=index_sort_c(i)
      ENDDO
!
!.....Get sorted global j_neigh_c_tmp to get nji
!
      DO ne=1,nelem
         j0=i_neigh_c_tmp(ne)-1
         DO j=i_neigh_c_tmp(ne),i_neigh_c_tmp(ne+1)-1
            k=index_sort0(j)
            temp_sort(j-j0)=j_neigh_c_tmp(j0+k)
         ENDDO
         DO j=i_neigh_c_tmp(ne),i_neigh_c_tmp(ne+1)-1
            j_neigh_c_tmp(j)=temp_sort(j-j0)
         ENDDO
         DO j=i_neigh_c_tmp(ne),i_neigh_c_tmp(ne+1)-1
            k=index_sort0(j)
            temp_sort(j-j0)=j_nbcon_c_tmp(j0+k)
         ENDDO
         DO j=i_neigh_c_tmp(ne),i_neigh_c_tmp(ne+1)-1
            j_nbcon_c_tmp(j)=temp_sort(j-j0)
         ENDDO
      ENDDO
!
      ALLOCATE(nji_c(maxmt_ncond))
      DO i=1,ncell_cond
         ne=i
         j0=i_neigh_c_tmp(ne)-1
         j3=i_neigh_c(i)-1
         DO j=i_neigh_c(i),i_neigh_c(i+1)-1
            IF(nbcon_c(j).eq.0)THEN
               k=j_neigh_c_tmp(j-j3+j0)
               j2=i_neigh_c_tmp(k)-1
               DO j1=i_neigh_c_tmp(k),i_neigh_c_tmp(k+1)-1
                  IF (j_neigh_c_tmp(j1).eq.ne) then
                      nji_c(j)=j1-j2
                      GOTO 200
                  ENDIF
               ENDDO
200            CONTINUE
            ELSE
!              nji_c(j,i)=0
            ENDIF
         ENDDO
      ENDDO
!
      maxmt_ps=i_neigh_c(ncell_ps+1)-1
!
      END SUBROUTINE subdomain_info_solid_ser
!
