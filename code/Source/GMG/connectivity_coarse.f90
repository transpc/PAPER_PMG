      Subroutine connectivity_coarse_nnz(nnode,nnode1,nnz,nnzi,ia,ja,iai,jai,  &
                 iar,jar,nnz1,jmax)

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * !
! find the connectivity - neibor_c - of coarse grid                             ! 
! for each node I, fine all nodes J such that                                   !
! Ac(I,J) = R(I,k)*A(k,l)*P(l,J) =/ 0 , A(k,l) on fine -grid.                   !
! because the grid is non-nested => J is a litle more than all neibor of I      !
!   desiged by Sang-Ha - May-2018                                               !
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
!  for each node I :                                                            !
!  find all nodes -k on fine-grid st: R(I,k)=/0 -> k is neighbor node of I      !
!  find nodes l on fine-grid st: A(k,l) =/ 0                                    !
!  find node J on coarse gird st: P(l,J) =/ 0                                   !
!  all matrix: R(restriction),A ,P(interpolation) stored in CSR format          !
! input:                                                                        !
! matrix A, matrix R, matrix P                                                  !
! output:                                                                       !
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! 
      implicit None
      
! input
      integer(4) nnode,nnode1,nnz,nnzi
      integer(4) ia(nnode+1),ja(nnz)
      integer(4) iai(nnode+1),jai(nnzi)
      integer(4) iar(nnode1+1),jar(nnzi)

! output
      integer nnz1,jmax

! free 
      integer i,j,k,l,i1,i2,k1,k2,l1,l2,j1,j2
      integer nnd,id,id1,id2
      integer,dimension(:),allocatable:: imark
    
! ----
      allocate(imark(nnode1))
      imark = 0
      nnz1 = 0
      jmax = 0
      
      DO i = 1,nnode1
          
        nnd = 0  
! ---for node i   
        
         i1 = iar(i)
         i2 = iar(i+1)-1		 
         do id = i1,i2
             k = jar(id)              ! R(I,k)
             k1 = ia(k)
             k2 = ia(k+1)-1
             do id1 = k1,k2
                 l = ja(id1)          ! A(k,l)
                 j1 = iai(l)
                 j2 = iai(l+1)-1
                 do id2 = j1,j2
                     j = jai(id2)     ! P(l,J)
                     
                     if(imark(j).ne.i) then
                         nnd = nnd + 1 
!                         ni(nnd) = j
                         imark(j) = i
                     end if
                     
                 end do       ! id2 loop
                 
             end do           ! id1 loop
             
         end do               ! id loop
         
         jmax = MAX(jmax,nnd)
        nnz1 = nnz1 + nnd   
        
      END DO
      
      Deallocate(imark)   
             	
!      write(*,*)'nnz on coarse grid = ',nnz1
! --------------------------------------
      
      return
      
    End	  
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
    
    Subroutine connectivity_coarse_CSR(jmax,nnode,nnode1,nnz,nnzi,ia,ja,iai,jai,  &
                 iar,jar,nnz1,ia1,ja1,ju1)

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * !
! find the connectivity - neibor_c - of coarse grid                             ! 
! for each node I, fine all nodes J such that                                   !
! Ac(I,J) = R(I,k)*A(k,l)*P(l,J) =/ 0 , A(k,l) on fine -grid.                   !
! because the grid is non-nested => J is a litle more than all neibor of I      !
!   desiged by Sang-Ha - May-2018                                               !
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
!  for each node I :                                                            !
!  find all nodes -k on fine-grid st: R(I,k)=/0 -> k is neighbor node of I      !
!  find nodes l on fine-grid st: A(k,l) =/ 0                                    !
!  find node J on coarse gird st: P(l,J) =/ 0                                   !
!  all matrix: R(restriction),A ,P(interpolation) stored in CSR format          !
! input:                                                                        !
! matrix A, matrix R, matrix P                                                  !
! output:                                                                       !
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! 
      implicit None
      
! input
      integer jmax,nnz1
      integer nnode,nnode1,nnz,nnzi
      integer ia(nnode+1),ja(nnz)
      integer iai(nnode+1),jai(nnzi)
      integer iar(nnode1+1),jar(nnzi)

! output
      integer ia1(nnode1+1), ja1(nnz1), ju1(nnode1)

! free 
      integer i,j,k,l,i1,i2,k1,k2,l1,l2,j1,j2
      integer nnd,id,nn,id1,id2
      integer ni(jmax)
      integer,dimension(:),allocatable:: imark
    
! ----
      allocate(imark(nnode1))
      imark = 0
      ia1(1) = 1
      
      DO i = 1,nnode1
          
        nnd = 0  
! ---for node i   
        
         i1 = iar(i)
         i2 = iar(i+1)-1		 
         do id = i1,i2
             k = jar(id)              ! R(I,k)
             k1 = ia(k)
             k2 = ia(k+1)-1
             do id1 = k1,k2
                 l = ja(id1)          ! A(k,l)
                 j1 = iai(l)
                 j2 = iai(l+1)-1
                 do id2 = j1,j2
                     j = jai(id2)     ! P(l,J)
                     
                     if(imark(j).ne.i) then
                         nnd = nnd + 1 
                         ni(nnd) = j
                         imark(j) = i
                     end if
                     
                 end do       ! id2 loop
                 
             end do           ! id1 loop
             
         end do               ! id loop

!   reordering 
        call bubble_sort(nnd,ni)
        
! CSR format 
! note that ni is including i (not only neighbor nodes)
        
        nn = ia1(i)
        ia1(i+1) = ia1(i) + nnd
        
        DO j = 1, nnd
            j1 = ni(j)
            j2 = nn+j-1
            ja1(j2) = j1
            IF(j1.EQ.i) ju1(i) = j2
        ENDDO
! - - - - - - -           
        
      END DO
      
      Deallocate(imark)  
      
      IF(ia1(nnode1+1).NE.(nnz1+1)) THEN
          WRITE(*,*)'PMG error in connectivity_coarse, ia1/nnz',ia1(nnode1+1),nnz1+1
          STOP
      ENDIF
      
             	
! --------------------------------------
      
      return
      
    End	  
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
    
		