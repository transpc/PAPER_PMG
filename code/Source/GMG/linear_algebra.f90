subroutine linalg_invM(n, a, option)
   !% ============================================
   !% Compute the inverse of a real matrix using
   !% LAPACK library
   !% input: - [n] : size of matrix a[n x n]
   !%          [option]: <1> - general matrix
   !%                    <2> - symmetric matrix
   !%                    <3> - symmetric (Hermitian)
   !%                          positive-definite matrix
   !% output: [a](n, n): overwritten
   !% ============================================
   implicit none
   !% use math kernel library
!   include 'mkl.fi'
   !% in/out
   integer(4), intent(in)    :: n, option
   real(8)   , intent(inout) :: a(n,*)
   !% working variables
   character(1) :: uplo
   integer(4)  :: m, lda, lwork, info
   integer(4), allocatable :: ipiv(:)
   real(8)   , allocatable :: work(:)
   
   !% initialize parameters/
   lwork = 16*n
   m     = n
   lda   = n
   allocate (ipiv(n), work(lwork))
   
   !% general matrix
   !% compute the LU factrization of a matrix /
   !% Syntax from Math kernel Libraries
   !< call dgetrf( m, n, a, lda, ipiv, info ) >
   call dgetrf_new(m, n, a, lda, ipiv, info)
   !% compute the inverse of a matrix /
   !% Syntax from Math kernel Libraries
   !< call dgetri( n, a, lda, ipiv, work, lwork, info ) >
   call dgetri_new(n, a, lda, ipiv, work, lwork, info)
   
   !% release allocated memory /
   deallocate (ipiv, work)
   return
    end subroutine
    
! = = = = = = = = = = = = = = = = = = = = = = = = = 

! inverse matrix by using GS elimination

subroutine matrix_inverse_GS(n,A)
      use omp_lib
      USE MD_OpenMP
  implicit none
  ! input
  integer :: n
  real(8) :: A(n,n)
  ! tmp
  real(8) ::  Ide(n,n)
  integer :: i, j, k
  real(8) :: pivot, temp
  ! Initialize the identity matrix I
  Ide = 0.0
  do i = 1, n
     Ide(i, i) = 1.0
  end do
  ! Perform Gauss-Jordan elimination
  do i = 1, n
     pivot = A(i, i)
     A(i, :) = A(i, :) / pivot
     Ide(i, :) = Ide(i, :) / pivot
    !$omp PARALLEL DO PRIVATE(temp)
     do j = 1, n
        if (j /= i) then
           temp = A(j, i)
           A(j, :) = A(j, :) - temp * A(i, :)
           Ide(j, :) = Ide(j, :) - temp * Ide(i, :)
        end if
     end do
    !$omp END PARALLEL DO
	
  end do
! out:
   A = Ide
  !
  return
end 

! == = = = = = = =  = = = = = = = = = = = = = = = = = = = = 
! inverse matrix by using GS elimination

subroutine matrix_inverse_GS_n(n,A)
      use omp_lib
      USE MD_OpenMP
  implicit none
  ! input
  integer :: n
  real(8) :: A(n,n)
  ! tmp
  real(8) ::  Ide(n,n)
  integer :: i, j, k
  real(8) :: pivot, temp
  ! Initialize the identity matrix I
  Ide = 0.0
  do i = 1, n
     Ide(i, i) = 1.0
  end do
  ! Perform Gauss-Jordan elimination
  do i = 1, n
     pivot = A(i, i)

!  !$omp parallel private(j,temp,k)
!     !$omp do
     do j = 1, n
        A(i, j) = A(i, j) / pivot
        Ide(i, j) = Ide(i, j) / pivot
     end do
!     !$omp end do

!     !$omp do
     do j = 1, n
        if (j /= i) then
           temp = A(j, i)
           do k = 1, n
              A(j, k) = A(j, k) - temp * A(i, k)
              Ide(j, k) = Ide(j, k) - temp * Ide(i, k)
           end do
        end if
     end do

!     !$omp end do
!  !$omp end parallel
  end do
  
! out:
   A = Ide
  !
  return
end 

