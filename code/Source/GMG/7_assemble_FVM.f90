   SUBROUTINE assemble_FVM(icase,nnz0c,source,au0c, diag)!, off_diag)
    
! reading matrix on finest mesh 
! it may not necessary for CUPID code 
!  

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * !
! assemble global matrix and RHS
! ---
   USE MD_parameter, ONLY: nf_max
   USE MD_geometry, ONLY: nelem,coord
   USE MD_matrix
   USE Ztimecon, ONLY : itim
   USE MD_MPI, ONLY: nintf,myrank
   USE MD_MG_matrix, ONLY: diagt

!----------------------------------!
   IMPLICIT NONE

   INTEGER(4) icase
   INTEGER ie,i1,i2,nd,j,id,i,nnz0c,k,j1,j2
   REAL(8) :: source(nintf), au0c (nnz0c), diag(nintf)
   REAL(8) xtmp,xtmp1,tmp,tmp2

!      
! --- stiffness matrix and RHS ---

  IF(icase.EQ.1) THEN

   b = 0.d0
   au = 0.d0
            
   DO ie = 1,nintf
!     
      b(ie) = source(ie)
      
      diagr(ie) = diag(ie)

      diagt(ie) = 1.d0/diag(ie)
          
   END DO
   
   do i = 1, nnz0c
      au(i) = au0c(i)
   end do
!   end do

!  
   ELSE

   b = 0.d0
   DO ie = 1,nintf
      b(ie) = source(ie)

   END DO

   ENDIF
!
 
    
   
   RETURN
   END
