      SUBROUTINE Prep_fine_P
! 
	  
! * * * * * * * * * * * * * * * * * * * * * * * *  * * * * * * * * * !

! ---
      USE MD_geometry,  ONLY: nnodegl,nnode,nelem, num_neigh,e_neigh, coord
      USE MD_matrix,    ONLY: nnz,ia,ja,ju,au,b,u,ut, diagr,alu,            &
	                &   nnz_l,nnz_u,ia_l,ia_u,ja_l,ja_u ,                   &
					&   iend,diag_l,diag_u,alu_l,alu_u 
      
      USE MD_parameter, ONLY: ndim,nf_max
      USE MD_MPI ,      ONLY: nnbd,nbdom,nnsend,nnrecv,spt,rpt,sintf,rintf,nintf,myrank
      USE MD_MPI_MG,    ONLY: siaf,riaf,nmaxgl
      USE MD_MG_matrix, ONLY: diagt
!
!-------------------
      IMPLICIT NONE

!!DEC$IF defined (mpi_flag)
!      INCLUDE 'mpif.h'
!!DEC$ENDIF
!!DEC$IF defined (mpi_flag)
!      INTEGER(4)::status(mpi_status_size),tag,ierr
!!DEC$ENDIF

      REAL(8) pi,tmp,tmp1
      INTEGER(4) i,j,l,nd,nnzss,nnzrr,i1,i2,j1,j2,ni(nf_max)
      INTEGER(4)::alstatus

! ---     
      ALLOCATE(ut(nelem),diagt(nelem))
      ut = 0.d0
      diagt = 0.d0
     
! --- 
! 
     ALLOCATE(diagr(nelem))
     IF(ALLOCATED(ia))     DEALLOCATE(ia)
     IF(ALLOCATED(ja))     DEALLOCATE(ja)
     IF(ALLOCATED(ju))     DEALLOCATE(ju)
     IF(ALLOCATED(au))     DEALLOCATE(au)     
      ALLOCATE(u(nnode),b(nnode))
! 
      diagr = 0.d0
      u = 0.d0
      b = 0.d0
!      
	  
! step 1: temporary setup for each processor. 
      nnz = sum(num_neigh(1:nelem))+nelem       ! FVM
! 	  
      allocate(au(nnz),ia(nelem+1),ja(nnz),ju(nelem))
      ia = 0
      ja = 0
      ju = 0
      au = 0.d0
	  
      call csr_FVM_P(nnodegl,nelem,nf_max,num_neigh,e_neigh,ia,ja,ju,nnz)   


! step 2: csr for external nodes:
      
      nnsend = spt(nnbd+1)-1
      nnrecv = rpt(nnbd+1)-1
      nmaxgl = 1
!
      IF(nnbd.EQ.0) THEN
       nnsend = 0
       nnrecv = 0
      ENDIF
!/
      IF(min(nnsend,nnrecv) .LT. 0) THEN
          write(*,*)' negative nnrecv',nnsend,nnrecv
      ENDIF
!/        
      allocate(siaf(nnsend+1),riaf(nnrecv+1),stat=alstatus)
!/
     IF (alstatus/=0) THEN
         WRITE(*,*)'not allocated, siaf',myrank
     ENDIF
!/
      siaf(1) = 1
      Do i=1,nnsend
          j=sintf(i)
          nd = ia(j+1)-ia(j)
          siaf(i+1)=siaf(i)+nd
          nmaxgl = MAX(nmaxgl,nd)
      End do
      
      nnzss = siaf(nnsend+1)-1
      
!     
      riaf(1) = 1
      Do i=1,nnrecv
          j=rintf(i)
          nd = ia(j+1)-ia(j)
          riaf(i+1)=riaf(i)+nd
          nmaxgl = MAX(nmaxgl,nd)
      End do      
      
! check:
      IF((riaf(nnrecv+1)-riaf(1)).NE.(ia(nnode+1)-ia(nintf+1))) THEN
          WRITE(*,*)'error in pre_fine_P, ria',myrank
          STOP
      ENDIF      
      
      nnzrr = riaf(nnrecv+1)-1
      
! step 2: adding the csr for external nodes:

    IF(nnbd.NE.0) THEN
      call send_receive_csr(ndim,nmaxgl,nnbd,nnode,nnodegl,spt,rpt,sintf,rintf,          &
              nnzss,nnzrr,nnsend,nnrecv,siaf,riaf,nbdom,coord,nnz,ia,ja)
 
    ENDIF
    
! ALU 
    
!          if(isth == 2) then
      ALLOCATE(alu(nnz))
      alu = 0.d0
!          endif
    
!100   CONTINUE
! 
! deallocate:
    
    DEALLOCATE(coord)
    DEALLOCATE(num_neigh, e_neigh)
    
! pbicg-ali 

    allocate(iend(nelem))
    
    nnz_l=0
    nnz_u=0
    DO i=1,nelem
        iend(i) = ia(i+1)-1
        l=ju(i)-ia(i)
        nnz_l = nnz_l + l
        l = iend(i)-ju(i)
        nnz_u = nnz_u + l
    ENDDO
            
    allocate(diag_l(nelem), diag_u(nelem))
    allocate(ia_l(nelem+1), ia_u(nelem+1)) 
    allocate(ja_l(nnz_l), ja_u(nnz_u))
    allocate(alu_l(nnz_l), alu_u(nnz_u))
    
    diag_l = 0.d0
    diag_u = 0.d0
    alu_l = 0.d0
    alu_u = 0.d0
    ia_l = 0
    ia_u = 0
    ja_l = 0
    ja_u = 0
!
	  
      RETURN
    END
      
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
      subroutine CSR_FVM_P(nelemgl,nelem,nnd,iwk,iwork,ia,ja,ju,nnz)
      implicit none
! ---
      integer(4) nnd,nelem,nnz,nelemgl
      integer(4) iwk(nelem),iwork(nnd,nelem)

! --- out
      integer(4) ia(nelem+1),ja(nnz),ju(nelem)
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
            IF(id.EQ.ie) ju(ie) = k
            k = k+1
        ENDDO
        
		 ia(ie+1) = k
         
      End do
      
      IF(ia(nelem+1).ne.(nnz+1)) THEN
          WRITE(999,*)'csr-FVM error'
          STOP
      ENDIF
      
! --
	  
      return
      End

! = = = = = = = = = = = = = = = = = = = 
!
