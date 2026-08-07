!
      SUBROUTINE csr_neigh(maxmt,ncell,ia,ja,ju,iEND,i_neigh,neigh,nbcon)
!
!     This routine defines a CSR format for neigh array, and the neigh CSR 
!     array is defined only ONCE.
!
      USE Zmpi , ONLY:  maxmt_fluid
!
      IMPLICIT NONE
!      
!     input
      INTEGER ncell,maxmt
      INTEGER i_neigh(ncell+1),neigh(maxmt_fluid),nbcon(maxmt_fluid)
!     output
      INTEGER ia(ncell+1),ja(maxmt),ju(ncell),iEND(ncell)
!     local variables
      INTEGER i,k
      INTEGER j,j0,j1
!
      DO i=1,ncell
         j0=i_neigh(i)-1
         j1=ia(i)-1
         k=1
         DO j=i_neigh(i),i_neigh(i+1)-1
            IF(nbcon(j).ne.0) exit
            IF(neigh(j).gt.i) exit
            k=k+1
         ENDDO
!
!     lower part
!
         DO j=1,k-1
            ja(j+j1)=neigh(j+j0)
         ENDDO
!
!     diagonal
!
            j=k
            ja(j+j1)=i
            ju(i)=j+j1
!
!     upper part
!
         DO j=k,i_neigh(i+1)-i_neigh(i)
            IF(nbcon(j+j0).ne.0) exit
            ja(j+j1+1)=neigh(j+j0)
         ENDDO
!
!.....Get iend array
!
         DO j=ia(i+1)-1,ia(i),-1
            IF(ja(j).le.ncell) THEN
               iend(i)=j
               exit
            ENDIF
         ENDDO
      ENDDO          
!
      END SUBROUTINE csr_neigh
!
      SUBROUTINE csr_neigh_c(maxmt,ncell,ia,ja,ju,iEND,i_neigh,neigh,nbcon)
!
!     This routine defines a CSR format for neigh array, and the neigh CSR 
!     array is defined only ONCE.
!
      USE Zmpi      , ONLY: maxmt_ncond
!
      IMPLICIT NONE
!      
!     input
      INTEGER ncell,maxmt
      INTEGER i_neigh(ncell+1),neigh(maxmt_ncond),nbcon(maxmt_ncond)
!     output
      INTEGER ia(ncell+1),ja(maxmt),ju(ncell),iEND(ncell)
!     local variables
      INTEGER i,k
      INTEGER j,j0,j1
!
      DO i=1,ncell
         j0=i_neigh(i)-1
         j1=ia(i)-1
         k=1
         DO j=i_neigh(i),i_neigh(i+1)-1
            IF(nbcon(j).ne.0) exit
            IF(neigh(j).gt.i) exit
            k=k+1
         ENDDO
!
!     lower part
!
         DO j=1,k-1
            ja(j+j1)=neigh(j+j0)
         ENDDO
!
!     diagonal
!
            j=k
            ja(j+j1)=i
            ju(i)=j+j1
!
!     upper part
!
         DO j=k,i_neigh(i+1)-i_neigh(i)
            IF(nbcon(j+j0).ne.0) exit
            ja(j+j1+1)=neigh(j+j0)
         ENDDO
!
!.....Get iend array
!
         DO j=ia(i+1)-1,ia(i),-1
            IF(ja(j).le.ncell) THEN
               iend(i)=j
               exit
            ENDIF
         ENDDO
      ENDDO          
!
      END SUBROUTINE csr_neigh_c
