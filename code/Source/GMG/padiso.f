         
!      call pardiso_solve(neq,ia(1),ja(1),au(1),x(1),arhsu(1))    

!c***********************************************************************          
      subroutine pardiso_solve(n,ia,ja,a,x,b)
!c========================================================================
      use omp_lib
!      include 'mkl_pardiso.f77'
      integer n
      integer ia(1),ja(1)
      real*8 a(1),b(1),x(1)

!c---- Internal solver memory pointer for 64-bit architectures
!c---- integer*8 pt(64)
      integer*8 pt(64)
!c---- all other variables
      integer maxfct,mnum,mtype,phase,nrhs,error,msglvl
      integer iparm(64)
!c---- local variables
      integer i,idum
      real*8 waltime1,waltime2,ddum,mem

!c---- fill all arrays containing matrix data.
      nrhs=1
      maxfct=1
      mnum=1

!c---- set up PARDISO control parameters
      do i=1,64
      iparm(i)=0
      enddo

!c---- no solver default
      iparm(1)=1
!c---- fill-in reordering from METIS
      iparm(2)=2
!c---- numbers of processors
!c      iparm(3)=omp_get_max_threads()
      iparm(3)= 1  !   ncores

!c---- no iterative-direct algorithm
      iparm(4)=0 !  31!  0
!c---- no user fill-in reducing permutation
      iparm(5)=0
!c---- =0 solution on the first n compoments of x
      iparm(6)=0
!c---- default logical fortran unit number for output
      iparm(7)=16
!c---- numbers of iterative refinement steps 
      iparm(8)=9
!c---- not in use
      iparm(9)=0
!c---- perturbe the pivot elements with 1E-13
      iparm(10)=13
!c---- use nonsymmetric permutation and scaling MPS
      iparm(11)=1!0  !   1
!c---- not in use
      iparm(12)=0
      iparm(13)=1 !  0 saha----sonsymmtric
!c---- Output: number of perturbed pivots
      iparm(14)=0
!c---- not in use
      iparm(15)=0
      iparm(16)=0
      iparm(17)=0
!c---- output: number of nonzeros in the factor LU
      iparm(18)=-1
!c---- output: Mflops for LU factorization
      iparm(19)=-1
!c---- output: numbers of CG Iterations
      iparm(20)=0
!c---- initialize error flag
      error=0
!c---- print statistical information
      msglvl=0
!c---- symmetric positive definite matrix
      mtype=11  !saha---note for nonsymmetric matric
!c      mtype=-2
!c---- initiliaze the internal solver memory pointer
      do i=1,64
      pt(i)=0
      enddo

!c---- reordering and symbolic Factorization
      phase=11
!      write(*,*) ' the number cores=',iparm(3)

      call pardiso(pt,maxfct,mnum,mtype,phase,n,a,ia,ja,idum,
     +             nrhs,iparm,msglvl,ddum,ddum,error)
!c      WRITE(*,*) ' Reordering completed...'
      IF(error.NE.0) THEN
      WRITE(*,*) ' (Factorization1) The following ERROR was detected:',
     +           error
      STOP
      ENDIF
!c---- factorization.
      phase=22

      call pardiso(pt,maxfct,mnum,mtype,phase,n,a,ia,ja,idum,
     +             nrhs,iparm,msglvl,ddum,ddum,error)

!c      WRITE(*,*) ' Factorization completed...'
      IF(error.NE.0) THEN
      WRITE(*,*) ' (Factorization2) The following ERROR was detected:',
     +           error
      STOP
      ENDIF

!c---- back substitution and iterative refinement
      iparm(8)=2
      phase=33

      call pardiso(pt,maxfct,mnum,mtype,phase,n,a,ia,ja,idum,
     +             nrhs,iparm,msglvl,b,x,error)

!c---- termination and release of memory
      phase=-1

      call pardiso(pt,maxfct,mnum,mtype,phase,n,ddum,idum,idum,
     +             idum,nrhs,iparm,msglvl,ddum,ddum,error)

      
! 
!      write(*,*)'padiso solver is finish'
      
      return
      end
