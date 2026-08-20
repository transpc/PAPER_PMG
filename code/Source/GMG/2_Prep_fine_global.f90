!
      SUBROUTINE read_mesh_FVM
      
      USE Zzone, ONLY : ncell_fluid_all
      USE Znode, ONLY : n_node
      USE Zcoord1, ONLY : xloc_tmp
      
! NOTES: - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! reading the finest mesh information                                  !
! it includs the number of cell: nelem, total node: nnode              !
! the cell-connectivity, the node on each element (inode)              !
! and the coordinate of cell-center. (coord)                           !
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
!
      USE MD_parameter
      USE MD_geometry
      USE Zparam, ONLY : ndim_cupid => ndim
      USE Znode, ONLY : nmax_vertex
!
      IMPLICIT NONE
!
      INTEGER :: i,k,j
      INTEGER :: m, n, l
!
! ----
!      ndim = ndim_cupid
      nv_max = nmax_vertex
      nnode_mg = n_node
      nnode = nelem      ! notes this one

      ALLOCATE(coord(ndim,nelem))
      
      DO i = 1,nelem

         IF(ndim_cupid == 2) THEN
            coord(1:2,i) = xloc_tmp(i,1:2)
         ELSE
            coord(1:3,i) = xloc_tmp(i,1:3)
         ENDIF
         
      ENDDO
 
!
!
      RETURN
      END	  
    ! = = = = = = = = = = = = =  = = = = = = = = = = = = = = = = = = = = = !
    
    
      subroutine Prep_fine_FVM
	  
! * * * * * * * * * * * * * * * * * * * * * * * *  * * * * * * * * * !

! ---
      use MD_geometry, ONLY: nelem,neigh_mg,num_neigh_mg
      use MD_matrix, ONLY: nnz, ia,ja,ju
      use MD_parameter, only: nf_max
!
!-------------------
      implicit none
      INTEGER(4) alstatus,status
      CHARACTER(100) :: command

! ---     

!
      nnz = sum(num_neigh_mg(1:nelem))+nelem       ! FVM
! 	  
      allocate(ia(nelem+1),ja(nnz))
      ia = 0
      ja = 0
!      ju = 0
! notes that there is no JU.
      call csr_FVM(nelem,nf_max,num_neigh_mg,neigh_mg,ia,ja,nnz)

!      DEALLOCATE(ju)
!
!/ create folder "MG_tmp" here 
!/ for WINDOW: 
!      command = 'if not exist PMG_pre mkdir PMG_pre'
 !     call execute_command_line(command,status)
      call system('mkdir MG_tmp')
!      IF(status /= 0) THEN
!          STOP
!      ENDIF
!/
      return
    End
    
! = = = = = = = = = = 
      subroutine CSR_FVM(nelem,nnd,iwk,iwork,ia,ja,nnz)
      implicit none
! ---
      integer(4) nnd,nelem,nnz
      integer(4) iwk(nelem),iwork(nnd,nelem)

! --- out
      integer(4) ia(nelem+1),ja(nnz)
! ---
      integer(4) i,j,k,id,imax,itemp,ie,nd
      integer(4) ni(nnd+1)

! --------------------------------------------------------!
 
      ni = 0

! ---
          ia(1) = 1
      Do ie = 1,nelem
       nd = iwk(ie)
       ni(1:nd) = iwork(1:nd,ie)
       ni(nd+1) = ie
       
       nd = nd+1
       call bubble_sort(nd,ni)    
       
        k = ia(ie)
        DO j = 1, nd
            id = ni(j)
            ja(k) = id
!            IF(id.EQ.ie) ju(ie) = k
            k = k+1
        ENDDO
        
		 ia(ie+1) = k
         
      End do
      
      IF(ia(nelem+1).ne.(nnz+1)) THEN
          WRITE(*,*)'PMG error: csr-FVM (2_Prep_fine_global)'
          STOP
      ENDIF
	  
      return
    End
    
    

      

