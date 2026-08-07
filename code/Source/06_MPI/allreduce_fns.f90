!
!     Those routines provide a user friendly interface for the collection of data in MPI mode
!
!  Real form
!       allgather_vec_r     gatherv giving count and vector displacement
!       allgatherv_r        gatherv giving izone=0,1,2,3 0=fluid,1=solid,2=fluid_core,3-fuel_rod
!       gatherv_r           gatherv for root 0 giving izone=0,1,2,3 0=fluid,1=solid,2=fluid_core,3-fuel_rod
!       scatterv_r          scatterv 1D variable giving izone=0,1 information
!       scatterv_ndim_r     scatterv ndim variable giving izone=0,1 information
!       scatterv_csr_r      scatterv variable  giving izone=0,1 information
!       scatterv_node_csr_r scatterv CSR node variable 
!       allreduce_r         allreduce sum vector
!       allreducei_r        allreduce sum vector input=output
!       allreduce_max_r     allreduce max vector
!       reduce_max_r        reduce max vector for root 0
!       allreducei_max_r    allreduce max vector input=output
!       allreduce_min_r     allreduce min vector
!       allreducei_min_r    allreduce min vector input=output
!  Integer form
!       allgather_vec_i
!       allgather_i         allgather 1 element put in vector size np
!       allgatherv_i
!       gatherv_i
!       scatterv_csr_i
!       scatterv_node_csr_i
!       scatterv_i
!       allgatherv_csr_i
!       allreduce_i
!       allreducei_i
!       allreduce_max_i
!       allreducei_max_i
!       allreduce_l
!
!       bcast_r             broadcast 1 variable from rank 0
!       bcast_i             broadcast 1 variable from rank 0
!       bcast_l             broadcast 1 variable from rank 0
!       broadcast_i         broadcast 1 variable from rank
!
!       finalize_mpi       finalize mpi  
!       barrier_mpi        barrier mpi 
!
!------------------------------------------------------------------------------
!
      SUBROUTINE gather_vec_r(a,n,b,n_all,recv_count,disp)
!
!DEC$IF defined (mpi_flag)
      USE Zcore        , ONLY: np
!      
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
!.....Input
      INTEGER :: n,n_all
      INTEGER :: recv_count(np),disp(np)
      REAL(8) :: a(n)
!.....Output
      REAL(8) :: b(n_all)
!.....Local variables
      INTEGER :: i
      INTEGER root,ierr
!
      IF(np.gt.1) THEN
         root=0
         CALL MPI_GATHERV(a,n              ,MPI_DOUBLE_PRECISION, &
                          b,recv_count,disp,MPI_DOUBLE_PRECISION, &
                          root,MPI_COMM_WORLD,ierr)
      ELSE
         DO i=1,n_all
            b(i)=a(i)
         ENDDO
      ENDIF
!
!DEC$ELSE
      IMPLICIT NONE
!
!.....Input
      INTEGER :: n,n_all
      INTEGER :: recv_count(1),disp(1)
      REAL(8) :: a(n)
!.....Output
      REAL(8) :: b(n_all)
!.....Local variables
      INTEGER :: i
!
      DO i=1,n_all
         b(i)=a(i)
      ENDDO
!DEC$ENDIF
      END SUBROUTINE gather_vec_r
!
!------------------------------------------------------------------------------
!
      SUBROUTINE allgather_vec_r(a,n,b,n_all,recv_count,disp)
!
!DEC$IF defined (mpi_flag)
      USE Zcore        , ONLY: np
!      
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
!.....Input
      INTEGER :: n,n_all
      INTEGER :: recv_count(np),disp(np)
      REAL(8) :: a(n)
!.....Output
      REAL(8) :: b(n_all)
!.....Local variables
      INTEGER :: i
      INTEGER :: ierr
!
      IF(np.gt.1) THEN
         CALL MPI_ALLGATHERV(a,n              ,MPI_DOUBLE_PRECISION, &
                             b,recv_count,disp,MPI_DOUBLE_PRECISION, &
                             MPI_COMM_WORLD,ierr)
      ELSE
         DO i=1,n_all
            b(i)=a(i)
         ENDDO
      ENDIF 
!
!DEC$ELSE
      IMPLICIT NONE
!
!.....Input
      INTEGER :: n,n_all
      INTEGER :: recv_count(1),disp(1)
      REAL(8) :: a(n)
!.....Output
      REAL(8) :: b(n_all)
!.....Local variables
      INTEGER :: i
!
      DO i=1,n_all
         b(i)=a(i)
      ENDDO
!DEC$ENDIF
      END SUBROUTINE allgather_vec_r
!
!------------------------------------------------------------------------------
!
      SUBROUTINE allgather_vec_i(a,n,b,n_all,recv_count,disp)
!
!DEC$IF defined (mpi_flag)
      USE Zcore        , ONLY: np
!      
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
!.....Input
      INTEGER :: n,n_all
      INTEGER :: recv_count(np),disp(np)
      INTEGER :: a(n)
!.....Output
      INTEGER :: b(n_all)
!.....Local variables
      INTEGER :: i
      INTEGER :: ierr
!
      IF(np.gt.1) THEN
         CALL MPI_ALLGATHERV(a,n              ,MPI_INTEGER, &
                             b,recv_count,disp,MPI_INTEGER, &
                             MPI_COMM_WORLD,ierr)
      ELSE
         DO i=1,n_all
            b(i)=a(i)
         ENDDO
      ENDIF 
!
!DEC$ELSE
      IMPLICIT NONE
!
!.....Input
      INTEGER :: n,n_all
      INTEGER :: recv_count(1),disp(1)
      INTEGER :: a(n)
!.....Output
      INTEGER :: b(n_all)
!.....Local variables
      INTEGER :: i
!
      DO i=1,n_all
         b(i)=a(i)
      ENDDO
!DEC$ENDIF
      END SUBROUTINE allgather_vec_i
!
!------------------------------------------------------------------------------
!
      SUBROUTINE gather_vec_i(a,n,b,n_all,recv_count,disp)
!
!DEC$IF defined (mpi_flag)
      USE Zcore        , ONLY: np
!      
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
!.....Input
      INTEGER :: n,n_all
      INTEGER :: recv_count(np),disp(np)
      INTEGER :: a(n)
!.....Output
      INTEGER :: b(n_all)
!.....Local variables
      INTEGER :: i
      INTEGER root,ierr
!
      IF(np.gt.1) THEN
         root=0
         CALL MPI_GATHERV(a,n              ,MPI_INTEGER, &
                          b,recv_count,disp,MPI_INTEGER, &
                          root,MPI_COMM_WORLD,ierr)
      ELSE
         DO i=1,n_all
            b(i)=a(i)
         ENDDO
      ENDIF 
!
!DEC$ELSE
      IMPLICIT NONE
!
!.....Input
      INTEGER :: n,n_all
      INTEGER :: recv_count(1),disp(1)
      INTEGER :: a(n)
!.....Output
      INTEGER :: b(n_all)
!.....Local variables
      INTEGER :: i
!
      DO i=1,n_all
         b(i)=a(i)
      ENDDO
!DEC$ENDIF
      END SUBROUTINE gather_vec_i
!
!------------------------------------------------------------------------------
!
      SUBROUTINE gatherv_r(a,n,b,n_all,izone)
!
!     izone=0 fluid
!     izone=1 solid
!     izone=2 fluid_core
!     izone=3 fuel_rod
!      
!DEC$IF defined (mpi_flag)
      USE IFCORE
      USE Zmpi      , ONLY: ncell_fluid1,ncell_fluid1_dsp,jjperm, &
                            ncell_fluid1_c,ncell_fluid1_dsp_c,jjperm_c
      USE Zcore     , ONLY: np,myrank
      USE Zrv_mpi   , ONLY: ncell_fluid1_core,ncell_fluid1_core_dsp,jjperm_fluid_core, &
                            ncell_fuel1_rod,ncell_fuel1_rod_dsp,jjperm_fuel_rod

!      
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
!.....Input
      INTEGER n,n_all,izone
      REAL(8) a(n)
!.....Output
      REAL(8) b(n_all)
!.....Local variables
      INTEGER i,j
      INTEGER root,ierr
!
      REAL(8), ALLOCATABLE::a1(:)
!
      IF(n_all.eq.0) THEN
         WRITE(*,*) 'gatherv_r should not be called with zero data count'
         WRITE(*,*) 'To display subroutine calling sequence recompile with option -g -traceback'
         CALL tracebackqq(user_exit_code=-1)
         CALL MPI_FINALIZE(ierr)
         STOP
      ENDIF
!
      IF(np.gt.1) THEN
         root=0
         IF(myrank.eq.root) THEN
            ALLOCATE (a1(n_all))
         ELSE
            ALLOCATE (a1(1))
         ENDIF
!
         IF(izone.eq.0)THEN
            CALL MPI_GATHERV(a ,n,                            MPI_DOUBLE_PRECISION, &
                             a1,ncell_fluid1,ncell_fluid1_dsp,MPI_DOUBLE_PRECISION, &
                             root,MPI_COMM_WORLD,ierr)
            IF(myrank.eq.root) THEN
               DO i=1,n_all
                   j=jjperm(i)
                   b(j)=a1(i)
               ENDDO
            ENDIF
         ELSEIF(izone.eq.1)THEN
            CALL MPI_GATHERV(a ,n,                                MPI_DOUBLE_PRECISION, &
                             a1,ncell_fluid1_c,ncell_fluid1_dsp_c,MPI_DOUBLE_PRECISION, &
                             root,MPI_COMM_WORLD,ierr)
            IF(myrank.eq.root) THEN
               DO i=1,n_all
                  j=jjperm_c(i)
                  b(j)=a1(i)
               ENDDO
            ENDIF
         ELSEIF(izone.eq.2)THEN
            CALL MPI_GATHERV(a ,n,                                      MPI_DOUBLE_PRECISION, &
                             a1,ncell_fluid1_core,ncell_fluid1_core_dsp,MPI_DOUBLE_PRECISION, &
                             root,MPI_COMM_WORLD,ierr)
            IF(myrank.eq.root) THEN
               DO i=1,n_all
                  j=jjperm_fluid_core(i)
                  b(j)=a1(i)
               ENDDO
            ENDIF
         ELSEIF(izone.eq.3)THEN
            CALL MPI_GATHERV(a ,n,                                  MPI_DOUBLE_PRECISION, &
                             a1,ncell_fuel1_rod,ncell_fuel1_rod_dsp,MPI_DOUBLE_PRECISION, &
                             root,MPI_COMM_WORLD,ierr)
            IF(myrank.eq.root) THEN
               DO i=1,n_all
                  j=jjperm_fuel_rod(i)
                  b(j)=a1(i)
               ENDDO
            ENDIF
         ENDIF
!
         DEALLOCATE(a1)
!
      ELSE
!
         DO i=1,n_all
            b(i)=a(i)
         ENDDO
      ENDIF
!
!DEC$ELSE
      IMPLICIT NONE
!.....Input
      INTEGER n,n_all,izone
      REAL(8) a(n)
!.....Output
      REAL(8) b(n_all)
!.....Local variables
      INTEGER i
!
      DO i=1,n_all
         b(i)=a(i)
      ENDDO
!
!DEC$ENDIF
!
      END SUBROUTINE gatherv_r
!
!------------------------------------------------------------------------------
!
      SUBROUTINE gatherv_r_2d(a,lda,b,na,nb,izone)
!
!     izone=0 fluid
!     izone=1 solid
!     izone=2 fluid_core
!     izone=3 fuel_rod
!
      USE Zparam     , ONLY: ndim
      USE Zrv_hts_2d , ONLY: nr_2d
!
!DEC$IF defined (mpi_flag)
      USE IFCORE
      USE Zmpi       , ONLY: ncell_fluid1,ncell_fluid1_dsp,jjperm, &
                             ncell_fluid1_2d,ncell_fluid1_2d_dsp
      USE Zcore      , ONLY: np,myrank
      USE Zrv_mpi    , ONLY: ncell_fuel1_rod,ncell_fuel1_rod_dsp,jjperm_fuel_rod,       &
                             ncell_fuel1_rod_2d,ncell_fuel1_rod_2d_dsp
!      
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
!.....Input
      INTEGER :: lda,na,nb,izone
      REAL(8) :: a(lda,*)
!.....Output
      REAL(8) :: b(nb,*)
!.....Local variables
      INTEGER :: i,j,l,jj,j0,ip
      INTEGER :: root,ierr
!
      REAL(8),DIMENSION(:),ALLOCATABLE :: a1
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: aa
!
      IF(nb.eq.0) THEN
         WRITE(*,*) 'gatherv_r_2d should not be called with zero data count'
         WRITE(*,*) 'To display subroutine calling sequence recompile with option -g -traceback'
         CALL tracebackqq(user_exit_code=-1)
         CALL MPI_FINALIZE(ierr)
         STOP
      ENDIF
!
      IF(np.gt.1) THEN
         root=0
!
         IF(izone.eq.0)THEN
            IF(myrank.eq.root) THEN
               ALLOCATE (a1(nb*ndim))
            ELSE
               ALLOCATE (a1(1))
            ENDIF
            IF(lda.eq.na) THEN
               CALL MPI_GATHERV(a ,na*ndim,                            MPI_DOUBLE_PRECISION, &
                                a1,ncell_fluid1_2d,ncell_fluid1_2d_dsp,MPI_DOUBLE_PRECISION, &
                                root,MPI_COMM_WORLD,ierr)
            ELSE
               ALLOCATE(aa(na,ndim))
               DO l=1,ndim
                  DO i=1,na
                     aa(i,l)=a(i,l)
                  ENDDO
               ENDDO
               CALL MPI_GATHERV(aa,na*ndim,                            MPI_DOUBLE_PRECISION, &
                                a1,ncell_fluid1_2d,ncell_fluid1_2d_dsp,MPI_DOUBLE_PRECISION, &
                                root,MPI_COMM_WORLD,ierr)
               DEALLOCATE(aa)
            ENDIF
            IF(myrank.eq.root) THEN
               DO ip=1,np
                  j0=ncell_fluid1_dsp(ip)
                  jj=1+ncell_fluid1_2d_dsp(ip)
                  DO l=1,ndim
                     DO i=1,ncell_fluid1(ip)
                        j=jjperm(i+j0)
                        b(j,l)=a1(jj)
                        jj=jj+1 
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF
!
         ELSEIF(izone.eq.3)THEN
!
            IF(myrank.eq.root) THEN
               ALLOCATE (a1(nb*nr_2d))
            ELSE
               ALLOCATE (a1(1))
            ENDIF
            IF(lda.eq.na) THEN
               CALL MPI_GATHERV(a ,na*nr_2d,                                 MPI_DOUBLE_PRECISION, &
                                a1,ncell_fuel1_rod_2d,ncell_fuel1_rod_2d_dsp,MPI_DOUBLE_PRECISION, &
                                root,MPI_COMM_WORLD,ierr)
            ELSE
               ALLOCATE(aa(na,nr_2d))
               DO l=1,nr_2d
                  DO i=1,na
                     aa(i,l)=a(i,l)
                  ENDDO
               ENDDO
               CALL MPI_GATHERV(aa,na*nr_2d,                                 MPI_DOUBLE_PRECISION, &
                                a1,ncell_fuel1_rod_2d,ncell_fuel1_rod_2d_dsp,MPI_DOUBLE_PRECISION, &
                                root,MPI_COMM_WORLD,ierr)
               DEALLOCATE(aa)
            ENDIF
            IF(myrank.eq.root) THEN
               DO ip=1,np
                  j0=ncell_fuel1_rod_dsp(ip)
                  jj=1+ncell_fuel1_rod_2d_dsp(ip)
                  DO l=1,nr_2d
                     DO i=1,ncell_fuel1_rod(ip)
                        j=jjperm_fuel_rod(i+j0)
                        b(j,l)=a1(jj)
                        jj=jj+1 
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF
         ENDIF
!
         DEALLOCATE(a1)
!
      ELSE
!
         IF(izone.eq.0)THEN
!
            DO l=1,ndim
               DO i=1,nb
                  b(i,l)=a(i,l)
               ENDDO
            ENDDO
!
         ELSEIF(izone.eq.3)THEN
!
            DO l=1,nr_2d
               DO i=1,nb
                  b(i,l)=a(i,l)
               ENDDO
            ENDDO
!
         ENDIF
!
      ENDIF
!DEC$ELSE
      IMPLICIT NONE
!
!.....Input
      INTEGER :: lda,na,nb,izone
      REAL(8) :: a(lda,*)
!.....Output
      REAL(8) :: b(nb,*)
!.....Local variables
      INTEGER :: i,l
!
      IF(izone.eq.0)THEN
!
         DO l=1,ndim
            DO i=1,nb
               b(i,l)=a(i,l)
            ENDDO
         ENDDO
!
      ELSEIF(izone.eq.3)THEN
!
         DO l=1,nr_2d
            DO i=1,nb
               b(i,l)=a(i,l)
            ENDDO
         ENDDO
!
      ENDIF
!
!DEC$ENDIF
!
      END SUBROUTINE gatherv_r_2d
!
!------------------------------------------------------------------------------
!
      SUBROUTINE allgatherv_r(a,b,n,n_all,izone)
!
!     izone=0 fluid
!     izone=1 solid
!     izone=2 fluid_core
!     izone=3 fuel_rod
!
!DEC$IF defined (mpi_flag)
      USE IFCORE
      USE Zmpi      , ONLY: ncell_fluid1,ncell_fluid1_dsp,jjperm, &
                            ncell_fluid1_c,ncell_fluid1_dsp_c,jjperm_c
      USE Zcore     , ONLY: np
      USE Zrv_mpi   , ONLY: ncell_fluid1_core,ncell_fluid1_core_dsp,jjperm_fluid_core, &
                            ncell_fuel1_rod,ncell_fuel1_rod_dsp,jjperm_fuel_rod

!      
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
!.....Input
      INTEGER n,n_all,izone
      REAL(8) a(n)
!.....Output
      REAL(8) b(n_all)
!.....Local variables
      INTEGER i,j
      INTEGER ierr
!
      REAL(8), ALLOCATABLE::a1(:)
!
      IF(n_all.eq.0) THEN
         WRITE(*,*) 'allgatherv_r should not be called with zero data count'
         WRITE(*,*) 'To display subroutine calling sequence recompile with option -g -traceback'
         CALL tracebackqq(user_exit_code=-1)
         CALL MPI_FINALIZE(ierr)
         STOP
      ENDIF
!
      IF(np.gt.1) THEN
         ALLOCATE (a1(n_all))
!
         IF(izone.eq.0)THEN
            CALL MPI_ALLGATHERV(a ,n,                            MPI_DOUBLE_PRECISION, &
                                a1,ncell_fluid1,ncell_fluid1_dsp,MPI_DOUBLE_PRECISION, &
                                MPI_COMM_WORLD,ierr)
            DO i=1,n_all
               j=jjperm(i)
               b(j)=a1(i)
            ENDDO
         ELSEIF(izone.eq.1)THEN
            CALL MPI_ALLGATHERV(a ,n,                                MPI_DOUBLE_PRECISION, &
                                a1,ncell_fluid1_c,ncell_fluid1_dsp_c,MPI_DOUBLE_PRECISION, &
                                MPI_COMM_WORLD,ierr)
            DO i=1,n_all
               j=jjperm_c(i)
               b(j)=a1(i)
            ENDDO
         ELSEIF(izone.eq.2)THEN
            CALL MPI_ALLGATHERV(a ,n,                                      MPI_DOUBLE_PRECISION, &
                                a1,ncell_fluid1_core,ncell_fluid1_core_dsp,MPI_DOUBLE_PRECISION, &
                                MPI_COMM_WORLD,ierr)
            DO i=1,n_all
               j=jjperm_fluid_core(i)
               b(j)=a1(i)
            ENDDO
         ELSEIF(izone.eq.3)THEN
            CALL MPI_ALLGATHERV(a ,n,                                  MPI_DOUBLE_PRECISION, &
                                a1,ncell_fuel1_rod,ncell_fuel1_rod_dsp,MPI_DOUBLE_PRECISION, &
                                MPI_COMM_WORLD,ierr)
            DO i=1,n_all
               j=jjperm_fuel_rod(i)
               b(j)=a1(i)
            ENDDO
         ENDIF
         DEALLOCATE(a1)
!
      ELSE
!
         DO i=1,n_all
            b(i)=a(i)
         ENDDO
      ENDIF
!
!DEC$ELSE
      IMPLICIT NONE
!.....Input
      INTEGER n,n_all,izone
      REAL(8) a(n)
!.....Output
      REAL(8) b(n_all)
!.....Local variables
      INTEGER i
!
      DO i=1,n_all
         b(i)=a(i)
      ENDDO
!
!DEC$ENDIF
!
      END SUBROUTINE allgatherv_r
!
!------------------------------------------------------------------------------
!
      SUBROUTINE allgatherv_r_2d(a,lda,b,na,nb,izone)
!
!     izone=0 fluid
!     izone=1 solid
!     izone=2 fluid_core
!     izone=3 fuel_rod
!
      USE Zparam     , ONLY: ndim
!DEC$IF defined (mpi_flag)
      USE IFCORE
      USE Zmpi      , ONLY: ncell_fluid1,ncell_fluid1_dsp,jjperm, &
                             ncell_fluid1_2d,ncell_fluid1_2d_dsp
      USE Zcore     , ONLY: np

!      
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
!.....Input
      INTEGER lda,na,nb,izone
      REAL(8) a(lda,*)
!.....Output
      REAL(8) b(nb,*)
!.....Local variables
      INTEGER i,j,l,jj,j0,ip
      INTEGER ierr
!
      REAL(8),DIMENSION(:),ALLOCATABLE :: a1
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: aa
!
      IF(nb.eq.0) THEN
         WRITE(*,*) 'allgatherv_r_2d should not be called with zero data count'
         WRITE(*,*) 'To display subroutine calling sequence recompile with option -g -traceback'
         CALL tracebackqq(user_exit_code=-1)
         CALL MPI_FINALIZE(ierr)
         STOP
      ENDIF
!
      IF(np.gt.1) THEN
         ALLOCATE (a1(nb*ndim))
         IF(izone.eq.0)THEN
            IF(lda.eq.na) THEN
               CALL MPI_ALLGATHERV(a ,na*ndim,                            MPI_DOUBLE_PRECISION, &
                                   a1,ncell_fluid1_2d,ncell_fluid1_2d_dsp,MPI_DOUBLE_PRECISION, &
                                   MPI_COMM_WORLD,ierr)
            ELSE
               ALLOCATE(aa(na,ndim))
               DO l=1,ndim
                  DO i=1,na
                     aa(i,l)=a(i,l)
                  ENDDO
               ENDDO
               CALL MPI_ALLGATHERV(aa,na*ndim,                            MPI_DOUBLE_PRECISION, &
                                   a1,ncell_fluid1_2d,ncell_fluid1_2d_dsp,MPI_DOUBLE_PRECISION, &
                                   MPI_COMM_WORLD,ierr)
            ENDIF
            DO ip=1,np
               j0=ncell_fluid1_dsp(ip)
               jj=1+ncell_fluid1_2d_dsp(ip)
               DO l=1,ndim
                  DO i=1,ncell_fluid1(ip)
                     j=jjperm(i+j0)
                     b(j,l)=a1(jj)
                     jj=jj+1
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
         DEALLOCATE(a1)
!
      ELSE
!
         IF(izone.eq.0)THEN
            DO l=1,ndim
               DO i=1,nb
                  b(i,l)=a(i,l)
               ENDDO
            ENDDO
         ENDIF
      ENDIF
!
!DEC$ELSE
      IMPLICIT NONE
!.....Input
      INTEGER :: lda,na,nb,izone
      REAL(8) :: a(lda,*)
!.....Output
      REAL(8) :: b(nb,*)
!.....Local variables
      INTEGER :: i,l
!
      IF(izone.eq.0)THEN
!
         DO l=1,ndim
            DO i=1,nb
               b(i,l)=a(i,l)
            ENDDO
         ENDDO
!
      ENDIF
!
!DEC$ENDIF
!
      END SUBROUTINE allgatherv_r_2d
!
!------------------------------------------------------------------------------
!
      SUBROUTINE allgather_i(a,b)
!
!DEC$IF defined (mpi_flag)
      USE Zcore        , ONLY: np
!      
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
!.....Input
      INTEGER a
!.....Output
      INTEGER b(np)
!.....Local variables
      INTEGER ierr
!
      IF(np.gt.1) THEN
         CALL MPI_ALLGATHER(a,1,MPI_INTEGER,     &
                            b,1,MPI_INTEGER,     &
                            MPI_COMM_WORLD,ierr)
      ELSE
         b(1)=a
      ENDIF
!
!DEC$ELSE
      IMPLICIT NONE
!
!.....Input
      INTEGER a
!.....Output
      INTEGER b
!
      b=a
!DEC$ENDIF
      END SUBROUTINE allgather_i
!
!------------------------------------------------------------------------------
!
      SUBROUTINE scatterv_r(b,a,n,n_all,izone)
!
!DEC$IF defined (mpi_flag)
      USE Zmpi   ,ONLY: ncell_fluid1,ncell_fluid1_dsp,     &
                        ncell_fluid1_c,ncell_fluid1_dsp_c, &
                        jjperm,jjperm_c
      USE Zcore  ,ONLY: np,myrank
!      
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
!     input
      INTEGER n,n_all,izone
      REAL(8) a(n_all)
!     output
      REAL(8) b(n)
!     local variable
      INTEGER i,ii,ip
      INTEGER ierr
!     local arrays
      REAL(8), ALLOCATABLE :: a1(:)
!
      IF(n_all.eq.0) THEN
         WRITE(*,*) 'scatterv_r should not be called with zero data count'
         CALL MPI_FINALIZE(ierr)
         STOP
      ENDIF
!
      IF(np.gt.1) THEN
         ip=myrank+1
         IF(myrank.eq.0) THEN
            ALLOCATE (a1(n_all))
         ELSE
            ALLOCATE (a1(1))
         ENDIF
         IF(izone.eq.0)THEN
            IF(myrank.eq.0) THEN
               DO ii=1,n_all
                  i=jjperm(ii)
                  a1(ii)=a(i)
               ENDDO
            ENDIF
            CALL MPI_SCATTERV(a1,ncell_fluid1,ncell_fluid1_dsp,MPI_DOUBLE_PRECISION, &
                              b,ncell_fluid1(ip),              MPI_DOUBLE_PRECISION, &
                              0,MPI_COMM_WORLD,ierr)
         ELSEIF(izone.eq.1)THEN
            IF(myrank.eq.0) THEN
               DO ii=1,n_all
                  i=jjperm_c(ii)
                  a1(ii)=a(i)
               ENDDO
            ENDIF
            CALL MPI_SCATTERV(a1,ncell_fluid1_c,ncell_fluid1_dsp_c,MPI_DOUBLE_PRECISION, &
                              b ,ncell_fluid1_c(ip),               MPI_DOUBLE_PRECISION, &
                              0,MPI_COMM_WORLD,ierr)
         ENDIF
         DEALLOCATE(a1)
      ELSE
         DO i=1,n
            b(i)=a(i)
         ENDDO
      ENDIF
!
!DEC$ELSE
!      
      IMPLICIT NONE
!.....Input
      INTEGER n,n_all,izone
      REAL(8) a(n_all)
!.....Output
      REAL(8) b(n)
!.....Local variable
      INTEGER i
!
      DO i=1,n
         b(i)=a(i)
      ENDDO
!DEC$ENDIF
!
      END SUBROUTINE scatterv_r
!
!------------------------------------------------------------------------------
!
      SUBROUTINE scatterv_ndim_r(b,a,n,n_all,izone)
!
      USE Zparam ,ONLY: ndim
!DEC$IF defined (mpi_flag)
      USE Zmpi   ,ONLY: ncell_ndim_sz,ncell_ndim_dsp,      &
                        ncell_ndim_sz_c,ncell_ndim_dsp_c,  &
                        ncell_fluid1,ncell_fluid1_dsp,     &
                        ncell_fluid1_c,ncell_fluid1_dsp_c, &
                        jjperm,jjperm_c
      USE Zcore  ,ONLY: myrank,np
!      
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
!     input
      INTEGER n,n_all,izone
      REAL(8) a(n_all,ndim)
!     output
      REAL(8) b(n,ndim)
!     local variable
      INTEGER i,ii,ix,ip
      INTEGER jp,i0,j0,nn
      INTEGER ierr
!     local arrays
      REAL(8), ALLOCATABLE :: a1(:)
!
      IF(n_all.eq.0) THEN
         WRITE(*,*) 'scatterv_ndim_r should not be called with zero data count'
         CALL MPI_FINALIZE(ierr)
         STOP
      ENDIF
!
      IF(np.gt.1) THEN
      ip=myrank+1
      IF(izone.eq.0) THEN
         IF(myrank.eq.0) THEN
            ALLOCATE (a1(n_all*ndim))
            DO jp=1,np
               nn=ncell_fluid1(jp)
               j0=ncell_fluid1_dsp(jp)
               i0=ncell_ndim_dsp(jp)
               DO ix=1,ndim
                  DO ii=1,nn
                     i=jjperm(j0+ii)
                     i0=i0+1
                     a1(i0)=a(i,ix)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            ALLOCATE (a1(1))
         ENDIF
         CALL MPI_SCATTERV(a1,ncell_ndim_sz,ncell_ndim_dsp,MPI_DOUBLE_PRECISION, &
                           b ,ncell_ndim_sz(ip),           MPI_DOUBLE_PRECISION, &
                           0,MPI_COMM_WORLD,ierr)
         DEALLOCATE(a1)
      ELSEIF(izone.eq.1) THEN
         IF(myrank.eq.0) THEN
            ALLOCATE (a1(n_all*ndim))
            DO jp=1,np
               nn=ncell_fluid1_c(jp)
               j0=ncell_fluid1_dsp_c(jp)
               i0=ncell_ndim_dsp_c(jp)
               DO ix=1,ndim
                  DO ii=1,nn
                     i=jjperm_c(j0+ii)
                     i0=i0+1
                     a1(i0)=a(i,ix)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            ALLOCATE (a1(1))
         ENDIF
         CALL MPI_SCATTERV(a1,ncell_ndim_sz_c,ncell_ndim_dsp_c,MPI_DOUBLE_PRECISION, &
                           b ,ncell_ndim_sz_c(ip),             MPI_DOUBLE_PRECISION, &
                           0,MPI_COMM_WORLD,ierr)
         DEALLOCATE(a1)
      ENDIF
      ELSE
         DO ix=1,ndim
            DO i=1,n
               b(i,ix)=a(i,ix)
            ENDDO
         ENDDO
      ENDIF
!
!DEC$ELSE
!      
      IMPLICIT NONE
!.....Input
      INTEGER n,n_all,izone
      REAL(8) a(n_all,ndim)
!.....Output
      REAL(8) b(n,ndim)
!.....Local variable
      INTEGER i,ix
!
      DO ix=1,ndim
         DO i=1,n
            b(i,ix)=a(i,ix)
         ENDDO
      ENDDO
!DEC$ENDIF
!
      END SUBROUTINE scatterv_ndim_r
!
!------------------------------------------------------------------------------
!
      SUBROUTINE scatterv_ndim_fp_r(b,a,n,n_fp,n_all,izone)
!
      USE Zparam ,ONLY: ndim
!DEC$IF defined (mpi_flag)
      USE Zmpi   ,ONLY: ncell_ndim_sz,ncell_ndim_dsp,      &
                        ncell_ndim_sz_c,ncell_ndim_dsp_c,  &
                        ncell_fluid1,ncell_fluid1_dsp,     &
                        ncell_fluid1_c,ncell_fluid1_dsp_c, &
                        jjperm,jjperm_c
      USE Zcore  ,ONLY: myrank,np
!      
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
!     input
      INTEGER n,n_fp,n_all,izone
      REAL(8) a(n_all,ndim)
!     output
      REAL(8) b(n_fp,ndim)
!     local variable
      INTEGER i,ii,ix,ip
      INTEGER jp,i0,j0
      INTEGER ierr
!     local arrays
      REAL(8),DIMENSION(:),ALLOCATABLE :: a1,b1
!
      IF(n_all.eq.0) THEN
         WRITE(*,*) 'scatterv_ndim_r should not be called with zero data count'
         CALL MPI_FINALIZE(ierr)
         STOP
      ENDIF
!
      IF(np.gt.1) THEN
         ip=myrank+1
         IF(myrank.eq.0) THEN
            ALLOCATE (a1(n_all*ndim))
         ELSE
            ALLOCATE (a1(1))
         ENDIF
         ALLOCATE (b1(n*ndim))
!
         IF(izone.eq.0) THEN
            IF(myrank.eq.0) THEN
               DO jp=1,np
                  j0=ncell_fluid1_dsp(jp)
                  i0=ncell_ndim_dsp(jp)
                  DO ix=1,ndim
                     DO ii=1,ncell_fluid1(jp)
                        i=jjperm(j0+ii)
                        i0=i0+1
                        a1(i0)=a(i,ix)
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF
            CALL MPI_SCATTERV(a1,ncell_ndim_sz,ncell_ndim_dsp,MPI_DOUBLE_PRECISION, &
                              b1,ncell_ndim_sz(ip),           MPI_DOUBLE_PRECISION, &
                              0,MPI_COMM_WORLD,ierr)
         ELSEIF(izone.eq.1) THEN
            IF(myrank.eq.0) THEN
               DO jp=1,np
                  j0=ncell_fluid1_dsp_c(jp)
                  i0=ncell_ndim_dsp_c(jp)
                  DO ix=1,ndim
                     DO ii=1,ncell_fluid1_c(jp)
                        i=jjperm_c(j0+ii)
                        i0=i0+1
                        a1(i0)=a(i,ix)
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF
            CALL MPI_SCATTERV(a1,ncell_ndim_sz_c,ncell_ndim_dsp_c,MPI_DOUBLE_PRECISION, &
                              b1,ncell_ndim_sz_c(ip),             MPI_DOUBLE_PRECISION, &
                              0,MPI_COMM_WORLD,ierr)
         ENDIF
         j0=0
         DO ix=1,ndim
            DO i=1,n
               j0=j0+1
               b(i,ix)=b1(j0)
            ENDDO
         ENDDO
         DEALLOCATE(a1,b1)
      ELSE
         DO ix=1,ndim
            DO i=1,n
               b(i,ix)=a(i,ix)
            ENDDO
         ENDDO
      ENDIF
!
!DEC$ELSE
!      
      IMPLICIT NONE
!.....Input
      INTEGER n,n_fp,n_all,izone
      REAL(8) a(n_all,ndim)
!.....Output
      REAL(8) b(n_fp,ndim)
!.....Local variable
      INTEGER i,ix
!
      DO ix=1,ndim
         DO i=1,n
            b(i,ix)=a(i,ix)
         ENDDO
      ENDDO
!DEC$ENDIF
!
      END SUBROUTINE scatterv_ndim_fp_r
!
!------------------------------------------------------------------------------
!
      SUBROUTINE scatterv_csr_r(b,nb_l,a,nb_all,n_all,nb,i_neigh_csr,izone)
!
!DEC$IF defined (mpi_flag)
      USE Zmpi   ,ONLY: ncell_csr_sz,ncell_csr_dsp,     &
                        ncell_csr_sz_c,ncell_csr_dsp_c, &
                        jjperm,jjperm_c
      USE Zcore  ,ONLY: np,myrank
!      
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
!     input
      INTEGER nb_l,nb_all
      INTEGER n_all,nb,izone
      INTEGER i_neigh_csr(n_all+1)
      REAL(8) a(nb_all)
!     output
      REAL(8) b(nb_l)
!     local variable
      INTEGER i,j,j1,ii,ip
      INTEGER ierr
!     local arrays
      REAL(8), ALLOCATABLE :: a1(:)
!
      IF(nb.eq.0) THEN
         WRITE(*,*) 'scatterv_csr_r should not be called with zero data count'
         CALL MPI_FINALIZE(ierr)
         STOP
      ENDIF
!
      IF(np.gt.1) THEN
         ip=myrank+1
         IF(myrank.eq.0) THEN
            ALLOCATE (a1(nb))
         ELSE
            ALLOCATE (a1(1))
         ENDIF
         IF(izone.eq.0)THEN
            IF(myrank.eq.0) THEN
               j1=0
               DO ii=1,n_all
                  i=jjperm(ii)
                  DO j=i_neigh_csr(i),i_neigh_csr(i+1)-1
                     j1=j1+1
                     a1(j1)=a(j)
                  ENDDO
               ENDDO
            ENDIF
            CALL MPI_SCATTERV(a1,ncell_csr_sz,ncell_csr_dsp,MPI_DOUBLE_PRECISION, &
                              b ,ncell_csr_sz(ip),          MPI_DOUBLE_PRECISION, &
                              0,MPI_COMM_WORLD,ierr)
         ELSEIF(izone.eq.1) THEN
            IF(myrank.eq.0) THEN
               j1=0
               DO ii=1,n_all
                  i=jjperm_c(ii)
                  DO j=i_neigh_csr(i),i_neigh_csr(i+1)-1
                     j1=j1+1
                     a1(j1)=a(j)
                  ENDDO
               ENDDO
            ENDIF
            CALL MPI_SCATTERV(a1,ncell_csr_sz_c,ncell_csr_dsp_c,MPI_DOUBLE_PRECISION, &
                              b ,ncell_csr_sz_c(ip),            MPI_DOUBLE_PRECISION, &
                              0,MPI_COMM_WORLD,ierr)
         ENDIF
         DEALLOCATE(a1)
      ELSE
         DO i=1,n_all
            DO j=i_neigh_csr(i),i_neigh_csr(i+1)-1
               b(j)=a(j)
            ENDDO
         ENDDO
      ENDIF
!
!DEC$ELSE
!      
      IMPLICIT NONE
!.....Input
      INTEGER nb_l,nb_all
      INTEGER n_all,nb,izone
      INTEGER i_neigh_csr(n_all+1)
      REAL(8) a(nb_all)
!.....Output
      REAL(8) b(nb_l)
!.....Local variable
      INTEGER i,j
!
      DO i=1,n_all
         DO j=i_neigh_csr(i),i_neigh_csr(i+1)-1
            b(j)=a(j)
         ENDDO
      ENDDO
!DEC$ENDIF
!
      END SUBROUTINE scatterv_csr_r
!
!------------------------------------------------------------------------------
!
      SUBROUTINE scatterv_csr_i(b,nb_l,a,nb_all,n_all,nb,i_neigh_csr,izone)
!
!DEC$IF defined (mpi_flag)
      USE Zmpi   ,ONLY: ncell_csr_sz,ncell_csr_dsp,     &
                        ncell_csr_sz_c,ncell_csr_dsp_c, &
                        jjperm,jjperm_c
      USE Zcore  ,ONLY: np,myrank
!      
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
!.....Input
      INTEGER nb_l,nb_all
      INTEGER n_all,nb,izone
      INTEGER i_neigh_csr(n_all+1)
      INTEGER a(nb_all)
!.....Output
      INTEGER b(nb_l)
!.....Local variable
      INTEGER i,j,j1,ii,ip
      INTEGER ierr
!     local arrays
      INTEGER, ALLOCATABLE :: a1(:)
!
      IF(nb.eq.0) THEN
         WRITE(*,*) 'scatterv_csr_i should not be called with zero data count'
         CALL MPI_FINALIZE(ierr)
         STOP
      ENDIF
!
      IF(np.gt.1) THEN
         ip=myrank+1
         IF(izone.eq.0)THEN
            IF(myrank.eq.0) THEN
               ALLOCATE (a1(nb))
               j1=0
               DO ii=1,n_all
                  i=jjperm(ii)
                  DO j=i_neigh_csr(i),i_neigh_csr(i+1)-1
                     j1=j1+1
                     a1(j1)=a(j)
                  ENDDO
               ENDDO
            ELSE
               ALLOCATE (a1(1))
            ENDIF
            CALL MPI_SCATTERV(a1,ncell_csr_sz,ncell_csr_dsp,MPI_INTEGER, &
                              b ,ncell_csr_sz(ip),          MPI_INTEGER, &
                              0,MPI_COMM_WORLD,ierr)
         ELSEIF(izone.eq.1) THEN
            IF(myrank.eq.0) THEN
               ALLOCATE (a1(nb))
               j1=0
               DO ii=1,n_all
                  i=jjperm_c(ii)
                  DO j=i_neigh_csr(i),i_neigh_csr(i+1)-1
                     j1=j1+1
                     a1(j1)=a(j)
                  ENDDO
               ENDDO
            ELSE
               ALLOCATE (a1(1))
            ENDIF
            CALL MPI_SCATTERV(a1,ncell_csr_sz_c,ncell_csr_dsp_c,MPI_INTEGER, &
                              b ,ncell_csr_sz_c(ip),            MPI_INTEGER, &
                              0,MPI_COMM_WORLD,ierr)
         ENDIF
         DEALLOCATE(a1)
      ELSE
         DO i=1,n_all
            DO j=i_neigh_csr(i),i_neigh_csr(i+1)-1
               b(j)=a(j)
            ENDDO
         ENDDO
      ENDIF
!
!DEC$ELSE
!      
      IMPLICIT NONE
!.....Input
      INTEGER nb_l,nb_all
      INTEGER n_all,nb,izone
      INTEGER i_neigh_csr(n_all+1)
      INTEGER a(nb_all)
!.....Output
      INTEGER b(nb_l)
!.....Local variable
      INTEGER i,j
!
      DO i=1,n_all
         DO j=i_neigh_csr(i),i_neigh_csr(i+1)-1
            b(j)=a(j)
         ENDDO
      ENDDO
!DEC$ENDIF
!
      END SUBROUTINE scatterv_csr_i
!
!------------------------------------------------------------------------------
!
      SUBROUTINE scatterv_csr_fp_i(b,nb_l,a,nb_all,n_all,nb,i_neigh_csr,izone,maxmt_fp,n,i_neigh)
!
!DEC$IF defined (mpi_flag)
      USE Zmpi   ,ONLY: ncell_csr_sz,ncell_csr_dsp,     &
                        ncell_csr_sz_c,ncell_csr_dsp_c, &
                        jjperm,jjperm_c
      USE Zcore  ,ONLY: np,myrank
!      
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
!.....Input
      INTEGER nb_l,nb_all
      INTEGER n_all,nb,izone,maxmt_fp,n
      INTEGER i_neigh(n+1)
      INTEGER i_neigh_csr(n_all+1)
      INTEGER a(nb_all)
!.....Output
      INTEGER b(maxmt_fp)
!.....Local variable
      INTEGER i,j,j0,j1,ii,ip
      INTEGER ierr
!     local arrays
      INTEGER,DIMENSION(:),ALLOCATABLE :: a1,b1
!
      IF(nb.eq.0) THEN
         WRITE(*,*) 'scatterv_csr_i should not be called with zero data count'
         CALL MPI_FINALIZE(ierr)
         STOP
      ENDIF
!
      IF(np.gt.1) THEN
         ip=myrank+1
         IF(myrank.eq.0) THEN
            ALLOCATE (a1(nb))
         ELSE
            ALLOCATE (a1(1))
         ENDIF
         ALLOCATE (b1(nb_l))
         IF(izone.eq.0)THEN
            IF(myrank.eq.0) THEN
               j1=0
               DO ii=1,n_all
                  i=jjperm(ii)
                  DO j=i_neigh_csr(i),i_neigh_csr(i+1)-1
                     j1=j1+1
                     a1(j1)=a(j)
                  ENDDO
               ENDDO
            ENDIF
            CALL MPI_SCATTERV(a1,ncell_csr_sz,ncell_csr_dsp,MPI_INTEGER, &
                              b1,ncell_csr_sz(ip),          MPI_INTEGER, &
                              0,MPI_COMM_WORLD,ierr)
         ELSEIF(izone.eq.1) THEN
            IF(myrank.eq.0) THEN
               j1=0
               DO ii=1,n_all
                  i=jjperm_c(ii)
                  DO j=i_neigh_csr(i),i_neigh_csr(i+1)-1
                     j1=j1+1
                     a1(j1)=a(j)
                  ENDDO
               ENDDO
            ENDIF
            CALL MPI_SCATTERV(a1,ncell_csr_sz_c,ncell_csr_dsp_c,MPI_INTEGER, &
                              b1,ncell_csr_sz_c(ip),            MPI_INTEGER, &
                              0,MPI_COMM_WORLD,ierr)
         ENDIF
         j0=0
         DO i=1,n
            DO j=i_neigh(i),i_neigh(i+1)-1
               j0=j0+1
               b(j)=b1(j0)
            ENDDO
         ENDDO
         DEALLOCATE(a1,b1)
      ELSE
         DO i=1,n_all
            DO j=i_neigh_csr(i),i_neigh_csr(i+1)-1
               b(j)=a(j)
            ENDDO
         ENDDO
      ENDIF
!
!DEC$ELSE
!      
      IMPLICIT NONE
!.....Input
      INTEGER nb_l,nb_all
      INTEGER n_all,nb,izone,maxmt_fp,n
      INTEGER i_neigh(n+1)
      INTEGER i_neigh_csr(n_all+1)
      INTEGER a(nb_all)
!.....Output
      INTEGER b(nb_l)
!.....Local variable
      INTEGER i,j
!
      DO i=1,n_all
         DO j=i_neigh_csr(i),i_neigh_csr(i+1)-1
            b(j)=a(j)
         ENDDO
      ENDDO
!DEC$ENDIF
!
      END SUBROUTINE scatterv_csr_fp_i
!
!------------------------------------------------------------------------------
!
      SUBROUTINE scatterv_csr_i_nbcon0(b,nb_l,a,nb_all,n_all,nb,i_neigh_csr,izone)
!
!DEC$IF defined (mpi_flag)
      USE Zmpi   ,ONLY: ncell_csr_sz_nbcon0,ncell_csr_dsp_nbcon0,     &
                        ncell_csr_sz_c,ncell_csr_dsp_c, &
                        jjperm,jjperm_c
      USE Zcore  ,ONLY: myrank
!      
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
!.....Input
      INTEGER nb_l,nb_all
      INTEGER n_all,nb,izone
      INTEGER i_neigh_csr(n_all+1)
      INTEGER a(nb_all)
!.....Output
      INTEGER b(nb_l)
!.....Local variable
      INTEGER i,j,j1,ii,ip
      INTEGER ierr
!     local arrays
      INTEGER, ALLOCATABLE :: a1(:)
!
      IF(nb.eq.0) THEN
         WRITE(*,*) 'scatterv_csr_i should not be called with zero data count'
         CALL MPI_FINALIZE(ierr)
         STOP
      ENDIF
!
      ip=myrank+1
      IF(izone.eq.0)THEN
         IF(myrank.eq.0) THEN
            ALLOCATE (a1(nb))
            j1=0
            DO ii=1,n_all
               i=jjperm(ii)
               DO j=i_neigh_csr(i),i_neigh_csr(i+1)-1
                  j1=j1+1
                  a1(j1)=a(j)
               ENDDO
            ENDDO
         ELSE
            ALLOCATE (a1(1))
         ENDIF
         CALL MPI_SCATTERV(a1,ncell_csr_sz_nbcon0,ncell_csr_dsp_nbcon0,MPI_INTEGER, &
                           b ,ncell_csr_sz_nbcon0(ip),                 MPI_INTEGER, &
                           0,MPI_COMM_WORLD,ierr)
         DEALLOCATE(a1)
      ELSEIF(izone.eq.1) THEN
         IF(myrank.eq.0) THEN
            ALLOCATE (a1(nb))
            j1=0
            DO ii=1,n_all
               i=jjperm_c(ii)
               DO j=i_neigh_csr(i),i_neigh_csr(i+1)-1
                  j1=j1+1
                  a1(j1)=a(j)
               ENDDO
            ENDDO
         ELSE
            ALLOCATE (a1(1))
         ENDIF
         CALL MPI_SCATTERV(a1,ncell_csr_sz_c,ncell_csr_dsp_c,MPI_INTEGER, &
                           b ,ncell_csr_sz_c(ip),            MPI_INTEGER, &
                           0,MPI_COMM_WORLD,ierr)
         DEALLOCATE(a1)
      ENDIF
!DEC$ENDIF
!
      END SUBROUTINE scatterv_csr_i_nbcon0
!
!------------------------------------------------------------------------------
!
      SUBROUTINE scatterv_node_csr_i(b,nb_l,a,nb_all,n_all,nb,i_cell_node_csr)
!
!DEC$IF defined (mpi_flag)
      USE Zmpi   ,ONLY: ncell_node_csr_sz,ncell_node_csr_dsp,     &
                        jjperm
      USE Zcore  ,ONLY: np,myrank
!      
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
!     input
      INTEGER nb_l,nb_all
      INTEGER n_all,nb
      INTEGER i_cell_node_csr(n_all+1)
      INTEGER a(nb_all)
!     output
      INTEGER b(nb_l)
!     local variable
      INTEGER i,j,j1,ii,ip
      INTEGER ierr
!     local arrays
      INTEGER, ALLOCATABLE :: a1(:)
!
      IF(nb.eq.0) THEN
         WRITE(*,*) 'scatterv_node_csr_i should not be called with zero data count'
         CALL MPI_FINALIZE(ierr)
         STOP
      ENDIF
!
      IF(np.gt.1) THEN
         ip=myrank+1
         IF(myrank.eq.0) THEN
            ALLOCATE (a1(nb))
            j1=0
            DO ii=1,n_all
               i=jjperm(ii)
               DO j=i_cell_node_csr(i),i_cell_node_csr(i+1)-1
                  j1=j1+1
                  a1(j1)=a(j)
               ENDDO
            ENDDO
         ELSE
            ALLOCATE (a1(1))
         ENDIF
         CALL MPI_SCATTERV(a1,ncell_node_csr_sz,ncell_node_csr_dsp,MPI_INTEGER, &
                           b ,ncell_node_csr_sz(ip),               MPI_INTEGER, &
                           0,MPI_COMM_WORLD,ierr)
          DEALLOCATE(a1)
      ELSE
         DO i=1,n_all
            DO j=i_cell_node_csr(i),i_cell_node_csr(i+1)-1
               b(j)=a(j)
            ENDDO
         ENDDO
      ENDIF
!
!DEC$ELSE
!      
      IMPLICIT NONE
!.....Input
      INTEGER nb_l,nb_all
      INTEGER n_all,nb,izone
      INTEGER i_cell_node_csr(n_all+1)
      INTEGER a(nb_all)
!.....Output
      INTEGER b(nb_l)
!.....Local variable
      INTEGER i,j
!
      DO i=1,n_all
         DO j=i_cell_node_csr(i),i_cell_node_csr(i+1)-1
            b(j)=a(j)
         ENDDO
      ENDDO
!DEC$ENDIF
!
      END SUBROUTINE scatterv_node_csr_i
!
!------------------------------------------------------------------------------
!
      SUBROUTINE scatterv_node_csr_fp_i(b,nb_l,a,nb_all,n_all,nd_max,n_fp,nb,i_cell_node_csr,n,num_cell_node)
!
!DEC$IF defined (mpi_flag)
      USE Zmpi   ,ONLY: ncell_node_csr_sz,ncell_node_csr_dsp,     &
                        jjperm
      USE Zcore  ,ONLY: np,myrank
!      
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
!     input
      INTEGER nb_l,nb_all,nd_max,n_fp
      INTEGER n_all,nb,n
      INTEGER num_cell_node(n_fp)
      INTEGER i_cell_node_csr(n_all+1)
      INTEGER a(nb_all)
!     output
      INTEGER b(nd_max,n_fp)
!     local variable
      INTEGER i,j,j0,j1,ii,ip
      INTEGER ierr
!     local arrays
      REAL(8),DIMENSION(:),ALLOCATABLE :: a1,b1
!
      IF(nb.eq.0) THEN
         WRITE(*,*) 'scatterv_node_csr_fp_i should not be called with zero data count'
         CALL MPI_FINALIZE(ierr)
         STOP
      ENDIF
!
      IF(np.gt.1) THEN
         ip=myrank+1
         ALLOCATE (b1(nb_l))
         IF(myrank.eq.0) THEN
            ALLOCATE (a1(nb))
            j1=0
            DO ii=1,n_all
               i=jjperm(ii)
               DO j=i_cell_node_csr(i),i_cell_node_csr(i+1)-1
                  j1=j1+1
                  a1(j1)=a(j)
               ENDDO
            ENDDO
         ELSE
            ALLOCATE (a1(1))
         ENDIF
         CALL MPI_SCATTERV(a1,ncell_node_csr_sz,ncell_node_csr_dsp,MPI_INTEGER, &
                           b1,ncell_node_csr_sz(ip),               MPI_INTEGER, &
                           0,MPI_COMM_WORLD,ierr)
         j0=0
         DO i=1,n
            DO j=1,num_cell_node(i)
               j0=j0+1
               b(j,i)=b1(j0)
            ENDDO
         ENDDO
         DEALLOCATE(a1,b1)
      ELSE
         j0=0
         DO i=1,n_all
            DO j=1,num_cell_node(i)
               j0=j0+1
               b(j,i)=a(j0)
            ENDDO
         ENDDO
      ENDIF
!
!DEC$ELSE
!      
      IMPLICIT NONE
!.....Input
      INTEGER nb_l,nb_all,nd_max,n_fp
      INTEGER n_all,nb,n
      INTEGER num_cell_node(n_fp)
      INTEGER i_cell_node_csr(n_all+1)
      INTEGER a(nb_all)
!.....Output
      INTEGER b(nd_max,n_fp)
!.....Local variable
      INTEGER i,j,j0
!
      j0=0
      DO i=1,n_all
         DO j=1,num_cell_node(i)
            j0=j0+1
            b(j,i)=a(j0)
         ENDDO
      ENDDO
!DEC$ENDIF
!
      END SUBROUTINE scatterv_node_csr_fp_i
!
!------------------------------------------------------------------------------
!
      SUBROUTINE scatterv_node_csr_r(b,nb_l,a,nb_all,n_all,nb,i_cell_node_csr)
!
!DEC$IF defined (mpi_flag)
      USE Zmpi   ,ONLY: ncell_node_csr_sz,ncell_node_csr_dsp,     &
                        jjperm
      USE Zcore  ,ONLY: np,myrank
!      
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
!     input
      INTEGER nb_l,nb_all
      INTEGER n_all,nb
      INTEGER i_cell_node_csr(n_all+1)
      REAL(8) a(nb_all)
!     output
      REAL(8) b(nb_l)
!     local variable
      INTEGER i,j,j1,ii,ip
      INTEGER ierr
!     local arrays
      REAL(8),DIMENSION(:),ALLOCATABLE :: a1
!
      IF(nb.eq.0) THEN
         WRITE(*,*) 'scatterv_node_csr_r should not be called with zero data count'
         CALL MPI_FINALIZE(ierr)
         STOP
      ENDIF
!
      IF(np.gt.1) THEN
         ip=myrank+1
         IF(myrank.eq.0) THEN
            ALLOCATE (a1(nb))
            j1=0
            DO ii=1,n_all
               i=jjperm(ii)
               DO j=i_cell_node_csr(i),i_cell_node_csr(i+1)-1
                  j1=j1+1
                  a1(j1)=a(j)
               ENDDO
            ENDDO
         ELSE
            ALLOCATE (a1(1))
         ENDIF
         CALL MPI_SCATTERV(a1,ncell_node_csr_sz,ncell_node_csr_dsp,MPI_DOUBLE_PRECISION, &
                           b ,ncell_node_csr_sz(ip),               MPI_DOUBLE_PRECISION, &
                           0,MPI_COMM_WORLD,ierr)
         DEALLOCATE(a1)
      ELSE
         DO i=1,n_all
            DO j=i_cell_node_csr(i),i_cell_node_csr(i+1)-1
               b(j)=a(j)
            ENDDO
         ENDDO
      ENDIF
!
!DEC$ELSE
!      
      IMPLICIT NONE
!.....Input
      INTEGER nb_l,nb_all
      INTEGER n_all,nb
      INTEGER i_cell_node_csr(n_all+1)
      REAL(8) a(nb_all)
!.....Output
      REAL(8) b(nb_l)
!.....Local variable
      INTEGER i,j
!
      DO i=1,n_all
         DO j=i_cell_node_csr(i),i_cell_node_csr(i+1)-1
            b(j)=a(j)
         ENDDO
      ENDDO
!DEC$ENDIF
!
      END SUBROUTINE scatterv_node_csr_r
!
!------------------------------------------------------------------------------
!
      SUBROUTINE scatterv_node_csr_fp_r(b,nb_l,a,nb_all,n_all,nd_max,n_fp,nb,i_cell_node_csr,n,num_cell_node)
!
!DEC$IF defined (mpi_flag)
      USE Zmpi   ,ONLY: ncell_node_csr_sz,ncell_node_csr_dsp,     &
                        jjperm
      USE Zcore  ,ONLY: np,myrank
!      
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
!     input
      INTEGER nb_l,nb_all,nd_max,n_fp
      INTEGER n_all,nb,n
      INTEGER num_cell_node(n_fp)
      INTEGER i_cell_node_csr(n_all+1)
      REAL(8) a(nb_all)
!     output
      REAL(8) b(nd_max,n_fp)
!     local variable
      INTEGER i,j,j0,j1,ii,ip
      INTEGER ierr
!     local arrays
      REAL(8),DIMENSION(:),ALLOCATABLE :: a1,b1
!
      IF(nb.eq.0) THEN
         WRITE(*,*) 'scatterv_node_csr_fp_r should not be called with zero data count'
         CALL MPI_FINALIZE(ierr)
         STOP
      ENDIF
!
      IF(np.gt.1) THEN
         ip=myrank+1
         ALLOCATE (b1(nb_l))
         IF(myrank.eq.0) THEN
            ALLOCATE (a1(nb))
            j1=0
            DO ii=1,n_all
               i=jjperm(ii)
               DO j=i_cell_node_csr(i),i_cell_node_csr(i+1)-1
                  j1=j1+1
                  a1(j1)=a(j)
               ENDDO
            ENDDO
         ELSE
            ALLOCATE (a1(1))
         ENDIF
         CALL MPI_SCATTERV(a1,ncell_node_csr_sz,ncell_node_csr_dsp,MPI_DOUBLE_PRECISION, &
                           b1,ncell_node_csr_sz(ip),               MPI_DOUBLE_PRECISION, &
                           0,MPI_COMM_WORLD,ierr)
         j0=0
         DO i=1,n
            DO j=1,num_cell_node(i)
               j0=j0+1
               b(j,i)=b1(j0)
            ENDDO
         ENDDO
         DEALLOCATE(a1,b1)
      ELSE
         j0=0
         DO i=1,n_all
            DO j=1,num_cell_node(i)
               j0=j0+1
               b(j,i)=a(j0)
            ENDDO
         ENDDO
      ENDIF
!
!DEC$ELSE
!      
      IMPLICIT NONE
!.....Input
      INTEGER nb_l,nb_all,nd_max,n_fp
      INTEGER n_all,nb,n
      INTEGER num_cell_node(n_fp)
      INTEGER i_cell_node_csr(n_all+1)
      REAL(8) a(nb_all)
!.....Output
      REAL(8) b(nd_max,n_fp)
!.....Local variable
      INTEGER i,j,j0
!
      j0=0
      DO i=1,n_all
         DO j=1,num_cell_node(i)
            j0=j0+1
            b(j,i)=a(j0)
         ENDDO
      ENDDO
!DEC$ENDIF
!
      END SUBROUTINE scatterv_node_csr_fp_r
!
!------------------------------------------------------------------------------
!
      SUBROUTINE scatterv_i(b,a,n,n_all,izone)
!
!DEC$IF defined (mpi_flag)
      USE Zmpi   ,ONLY: ncell_fluid1,ncell_fluid1_dsp,     &
                        ncell_fluid1_c,ncell_fluid1_dsp_c, &
                        jjperm,jjperm_c
      USE Zcore  ,ONLY: np,myrank
!      
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
!     input
      INTEGER a(n_all)
!     output
      INTEGER b(n)
!     local variable
      INTEGER i,ii,n,n_all,izone,ip
      INTEGER ierr
!     local arrays
      INTEGER,ALLOCATABLE :: a1(:)
!
      IF(n_all.eq.0) THEN
         WRITE(*,*) 'scatterv_i should not be called with zero data count'
         CALL MPI_FINALIZE(ierr)
         STOP
      ENDIF
!
      IF(np.gt.1) THEN
      ip=myrank+1
      IF(izone.eq.0) THEN
         IF(myrank.eq.0) THEN
            ALLOCATE (a1(n_all))
            DO ii=1,n_all
               i=jjperm(ii)
               a1(ii)=a(i)
            ENDDO
         ELSE
            ALLOCATE (a1(1))
         ENDIF
         CALL MPI_SCATTERV(a1,ncell_fluid1,ncell_fluid1_dsp,MPI_INTEGER, &
                           b ,ncell_fluid1(ip),             MPI_INTEGER, &
                           0,MPI_COMM_WORLD,ierr)
         DEALLOCATE(a1)
      ELSEIF(izone.eq.1) THEN
         IF(myrank.eq.0) THEN
            ALLOCATE (a1(n_all))
            DO ii=1,n_all
               i=jjperm_c(ii)
               a1(ii)=a(i)
            ENDDO
         ELSE
            ALLOCATE (a1(1))
         ENDIF
         CALL MPI_SCATTERV(a1,ncell_fluid1_c,ncell_fluid1_dsp_c,MPI_INTEGER, &
                           b ,ncell_fluid1_c(ip),               MPI_INTEGER, &
                           0,MPI_COMM_WORLD,ierr)
         DEALLOCATE(a1)
      ENDIF
      ELSE
         DO i=1,n
            b(i)=a(i)
         ENDDO
      ENDIF
!
!DEC$ELSE
!      
      IMPLICIT NONE
!.....Input
      INTEGER n,n_all,izone
      INTEGER a(n_all)
!.....Output
      INTEGER b(n)
!.....Local variable
      INTEGER i
!
      DO i=1,n
         b(i)=a(i)
      ENDDO
!DEC$ENDIF
!
      END SUBROUTINE scatterv_i
!
!------------------------------------------------------------------------------
!
      SUBROUTINE gatherv_i(a,n,b,n_all,izone)
!
!
!DEC$IF defined (mpi_flag)
      USE IFCORE
      USE Zmpi      ,ONLY: ncell_fluid1,ncell_fluid1_dsp,jjperm,      &
                           ncell_fluid1_c,ncell_fluid1_dsp_c,jjperm_c
      USE Zcore     ,ONLY: np,myrank
      USE Zrv_mpi   ,ONLY: ncell_fluid1_core,ncell_fluid1_core_dsp,jjperm_fluid_core
      USE Zrv_mpi   , ONLY: ncell_fuel1_rod,ncell_fuel1_rod_dsp,jjperm_fuel_rod
!      
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
!.....Input
      INTEGER :: n,n_all,izone
      INTEGER :: a(n)
!.....Output
      INTEGER :: b(n_all)
!.....Local variables
      INTEGER :: i,ii
      INTEGER :: root,ierr
!
      INTEGER,DIMENSION(:),ALLOCATABLE :: a1
!
      IF(n_all.eq.0) THEN
         WRITE(*,*) 'gatherv_i should not be called with zero data count'
         WRITE(*,*) 'To display subroutine calling sequence recompile with option -g -traceback'
         CALL tracebackqq(user_exit_code=-1)
         CALL MPI_FINALIZE(ierr)
         STOP
      ENDIF
!
      IF(np.gt.1) THEN
      root=0
      ALLOCATE (a1(n_all))
!
      IF(izone.eq.0)THEN
         CALL MPI_GATHERV(a ,n,                            MPI_INTEGER, &
                          a1,ncell_fluid1,ncell_fluid1_dsp,MPI_INTEGER, &
                          root,MPI_COMM_WORLD,ierr)
         IF(myrank.eq.root) THEN
            DO ii=1,n_all
               i=jjperm(ii)
               b(i)=a1(ii)
            ENDDO
         ENDIF
      ELSEIF(izone.eq.1)THEN
         CALL MPI_GATHERV(a ,n,                                MPI_INTEGER, &
                          a1,ncell_fluid1_c,ncell_fluid1_dsp_c,MPI_INTEGER, &
                          root,MPI_COMM_WORLD,ierr)
         IF(myrank.eq.root) THEN
            DO ii=1,n_all
               i=jjperm_c(ii)
               b(i)=a1(ii)
            ENDDO
         ENDIF
      ELSEIF(izone.eq.2)THEN
         CALL MPI_GATHERV(a ,n,                                      MPI_INTEGER, &
                          a1,ncell_fluid1_core,ncell_fluid1_core_dsp,MPI_INTEGER, &
                          root,MPI_COMM_WORLD,ierr)
         IF(myrank.eq.root) THEN
            DO ii=1,n_all
               i=jjperm_fluid_core(ii)
               b(i)=a1(ii)
            ENDDO
         ENDIF
      ELSEIF(izone.eq.3)THEN
         CALL MPI_GATHERV(a ,n,                                  MPI_INTEGER, &
                             a1,ncell_fuel1_rod,ncell_fuel1_rod_dsp,MPI_INTEGER, &
                             root,MPI_COMM_WORLD,ierr)
         IF(myrank.eq.root) THEN
            DO ii=1,n_all
                i=jjperm_fuel_rod(ii)
                b(i)=a1(ii)
            ENDDO
         ENDIF
      ENDIF
!
      DEALLOCATE(a1)
!
      ELSE
!
         DO i=1,n_all
            b(i)=a(i)
         ENDDO
!
      ENDIF
!
!DEC$ELSE
!
!      
      IMPLICIT NONE
!
!.....Input
      INTEGER :: n,n_all,izone
      INTEGER :: a(n)
!.....Output
      INTEGER :: b(n_all)
!.....Local variables
      INTEGER :: i
!
      DO i=1,n_all
         b(i)=a(i)
      ENDDO
!DEC$ENDIF
!
      END SUBROUTINE gatherv_i
!
!------------------------------------------------------------------------------
!
      SUBROUTINE allgatherv_i(a,b,n,n_all,izone)
!
!DEC$IF defined (mpi_flag)
      USE IFCORE
      USE Zmpi      , ONLY: ncell_fluid1,ncell_fluid1_dsp,jjperm,      &
                            ncell_fluid1_c,ncell_fluid1_dsp_c,jjperm_c
      USE Zcore     , ONLY: np
      USE Zrv_mpi   , ONLY: ncell_fluid1_core,ncell_fluid1_core_dsp,jjperm_fluid_core
      USE Zrv_mpi   , ONLY: ncell_fuel1_rod,ncell_fuel1_rod_dsp,jjperm_fuel_rod
!      
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
!     input
      INTEGER n,n_all,izone
      INTEGER a(n)
!     output
      INTEGER b(n_all)
!     local variables
      INTEGER i,ii
      INTEGER ierr
!
      INTEGER, ALLOCATABLE :: a1(:)
!
      IF(n_all.eq.0) THEN
         WRITE(*,*) 'allgatherv_i should not be called with zero data count'
         WRITE(*,*) 'To display subroutine calling sequence recompile with option -g -traceback'
         CALL tracebackqq(user_exit_code=-1)
         CALL MPI_FINALIZE(ierr)
         STOP
      ENDIF
!
      IF(np.gt.1) THEN
      ALLOCATE (a1(n_all))
!
      IF(izone.eq.0)THEN
         CALL MPI_ALLGATHERV(a,n,                             MPI_INTEGER, &
                             a1,ncell_fluid1,ncell_fluid1_dsp,MPI_INTEGER, &
                             MPI_COMM_WORLD,ierr)
         DO ii=1,n_all
             i=jjperm(ii)
             b(i)=a1(ii)
         ENDDO
      ELSEIF(izone.eq.1)THEN
         CALL MPI_ALLGATHERV(a,n,                                 MPI_INTEGER, &
                             a1,ncell_fluid1_c,ncell_fluid1_dsp_c,MPI_INTEGER, &
                             MPI_COMM_WORLD,ierr)
         DO ii=1,n_all
             i=jjperm_c(ii)
             b(i)=a1(ii)
         ENDDO
      ELSEIF(izone.eq.2)THEN
         CALL MPI_ALLGATHERV(a,n,                                       MPI_INTEGER, &
                             a1,ncell_fluid1_core,ncell_fluid1_core_dsp,MPI_INTEGER, &
                             MPI_COMM_WORLD,ierr)
         DO ii=1,n_all
             i=jjperm_fluid_core(ii)
             b(i)=a1(ii)
         ENDDO
      ELSEIF(izone.eq.3)THEN
         CALL MPI_ALLGATHERV(a ,n,                                  MPI_INTEGER, &
                             a1,ncell_fuel1_rod,ncell_fuel1_rod_dsp,MPI_INTEGER, &
                             MPI_COMM_WORLD,ierr)
         DO ii=1,n_all
             i=jjperm_fuel_rod(ii)
             b(i)=a1(ii)
         ENDDO
      ENDIF
      DEALLOCATE(a1)
      ELSE
!
         DO i=1,n_all
            b(i)=a(i)
         ENDDO
      ENDIF
!
!DEC$ELSE
      IMPLICIT NONE
!.....Input
      INTEGER n,n_all,izone
      INTEGER a(n)
!.....Output
      INTEGER b(n_all)
!.....Local variables
      INTEGER i
!
      DO i=1,n_all
         b(i)=a(i)
      ENDDO
!
!DEC$ENDIF
!
      END SUBROUTINE allgatherv_i
!
!------------------------------------------------------------------------------
!
      SUBROUTINE gatherv_csr_i(a,nb_l,b,nb_all,n_all,nb,i_neigh_csr,izone)
!
!DEC$IF defined (mpi_flag)
      USE Zmpi   ,ONLY: ncell_csr_sz,ncell_csr_dsp,     &
                        ncell_csr_sz_c,ncell_csr_dsp_c, &
                        jjperm,jjperm_c
      USE Zcore  ,ONLY: myrank
!      
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
!     input
      INTEGER nb_l,nb_all
      INTEGER n_all,nb,izone
      INTEGER i_neigh_csr(n_all+1)
      INTEGER a(nb_l)
!     output
      INTEGER b(nb_all)
!     local variables
      INTEGER i,j,ii,ip,j1
      INTEGER root,ierr
!     local arrays
      INTEGER, ALLOCATABLE :: a1(:)
!
      IF(n_all.eq.0) THEN
         WRITE(*,*) 'gatherv_csr_i should not be called with zero data count'
         CALL MPI_FINALIZE(ierr)
         STOP
      ENDIF
      ALLOCATE (a1(nb))
!
      root=0
      ip=myrank+1
      IF(izone.eq.0)THEN
         CALL MPI_GATHERV(a,ncell_csr_sz(ip),           MPI_INTEGER, &
                          a1,ncell_csr_sz,ncell_csr_dsp,MPI_INTEGER, &
                          root,MPI_COMM_WORLD,ierr)
         IF(myrank.eq.root) THEN
         j1=0
         DO ii=1,n_all
            i=jjperm(ii)
            DO j=i_neigh_csr(i),i_neigh_csr(i+1)-1
               j1=j1+1
               b(j)=a1(j1)
            ENDDO
         ENDDO
         ENDIF
      ELSEIF(izone.eq.1) THEN
         CALL MPI_GATHERV(a,ncell_csr_sz_c(ip),             MPI_INTEGER, &
                          a1,ncell_csr_sz_c,ncell_csr_dsp_c,MPI_INTEGER, &
                          root,MPI_COMM_WORLD,ierr)
         IF(myrank.eq.root) THEN
         j1=0
         DO ii=1,n_all
            i=jjperm_c(ii)
            DO j=i_neigh_csr(i),i_neigh_csr(i+1)-1
               j1=j1+1
               b(j)=a1(j1)
            ENDDO
         ENDDO
         ENDIF
      ENDIF
      DEALLOCATE(a1)
!DEC$ENDIF
!
      END SUBROUTINE gatherv_csr_i
!
!------------------------------------------------------------------------------
!
      SUBROUTINE allgatherv_csr_i(a,nb_l,b,nb_all,n_all,nb,i_neigh_csr,izone)
!
!DEC$IF defined (mpi_flag)
      USE Zmpi   ,ONLY: ncell_csr_sz,ncell_csr_dsp,     &
                        ncell_csr_sz_c,ncell_csr_dsp_c, &
                        jjperm,jjperm_c
      USE Zcore  ,ONLY: myrank
!      
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
!     input
      INTEGER nb_l,nb_all
      INTEGER n_all,nb,izone
      INTEGER i_neigh_csr(n_all+1)
      INTEGER a(nb_l)
!     output
      INTEGER b(nb_all)
!     local variables
      INTEGER i,j,ii,ip,j1
      INTEGER ierr
!     local arrays
      INTEGER, ALLOCATABLE :: a1(:)
!
      IF(n_all.eq.0) THEN
         WRITE(*,*) 'allgatherv_csr_i should not be called with zero data count'
         CALL MPI_FINALIZE(ierr)
         STOP
      ENDIF
      ALLOCATE (a1(nb))
!
      ip=myrank+1
      IF(izone.eq.0)THEN
         CALL MPI_ALLGATHERV(a,ncell_csr_sz(ip),           MPI_INTEGER, &
                             a1,ncell_csr_sz,ncell_csr_dsp,MPI_INTEGER, &
                             MPI_COMM_WORLD,ierr)
         j1=0
         DO ii=1,n_all
            i=jjperm(ii)
            DO j=i_neigh_csr(i),i_neigh_csr(i+1)-1
               j1=j1+1
               b(j)=a1(j1)
            ENDDO
         ENDDO
      ELSEIF(izone.eq.1) THEN
         CALL MPI_ALLGATHERV(a,ncell_csr_sz_c(ip),             MPI_INTEGER, &
                             a1,ncell_csr_sz_c,ncell_csr_dsp_c,MPI_INTEGER, &
                             MPI_COMM_WORLD,ierr)
         j1=0
         DO ii=1,n_all
            i=jjperm_c(ii)
            DO j=i_neigh_csr(i),i_neigh_csr(i+1)-1
               j1=j1+1
               b(j)=a1(j1)
            ENDDO
         ENDDO
      ENDIF
      DEALLOCATE(a1)
!DEC$ENDIF
!
      END SUBROUTINE allgatherv_csr_i
!
!------------------------------------------------------------------------------
!
      SUBROUTINE reducei_r(a,n)
!
!DEC$IF defined (mpi_flag)
      USE Zcore        , ONLY: np
!
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
!.....Input
      INTEGER :: n
      REAL(8) :: a(n)
!.....Local variables
      INTEGER :: i
      INTEGER :: root,ierr
      REAL(8) :: b(n)
!
      IF(np.gt.1) THEN
         root=0
         CALL MPI_REDUCE(a,b,n,MPI_DOUBLE_PRECISION,MPI_SUM, &
                         root,MPI_COMM_WORLD,ierr)
         DO i=1,n
            a(i)=b(i)
         ENDDO
      ENDIF
!DEC$ELSE
      IMPLICIT NONE
!
!.....Input
      INTEGER :: n
      REAL(8) :: a(n)
!
!DEC$ENDIF
!
      END SUBROUTINE reducei_r
!
!------------------------------------------------------------------------------
!
      SUBROUTINE reduce_r1(a,b)
!
!DEC$IF defined (mpi_flag)
      USE Zcore        , ONLY: np
!
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
!.....Input
      REAL(8) a
!.....Output
      REAL(8) b
!.....Local variables
      INTEGER root,ierr
!
      IF(np.gt.1) THEN
         root=0
         CALL MPI_REDUCE(a,b,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                         root,MPI_COMM_WORLD,ierr)
      ELSE
         a=b
      ENDIF
!DEC$ELSE
      IMPLICIT NONE
!
!.....Input
      REAL(8) a
!.....Output
      REAL(8) b
!
      a=b
!DEC$ENDIF
!
      END SUBROUTINE reduce_r1
!
!------------------------------------------------------------------------------
!
      SUBROUTINE allreduce_r(a,b,n)
!
!DEC$IF defined (mpi_flag)
      USE Zcore        , ONLY: np
!
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
!.....Input
      INTEGER :: n
      REAL(8) :: a(n)
!.....Output
      REAL(8) :: b(n)
!.....Local variables
      INTEGER :: i
      INTEGER :: ierr
!
      IF(np.gt.1) THEN
         CALL MPI_ALLREDUCE(a,b,n,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      ELSE
         DO i=1,n
            b(i)=a(i)
         ENDDO
      ENDIF
!DEC$ELSE
      IMPLICIT NONE
!
!.....Input
      INTEGER :: n
      REAL(8) :: a(n)
!.....Output
      REAL(8) :: b(n)
!.....Local variables
      INTEGER :: i
!
      DO i=1,n
         b(i)=a(i)
      ENDDO
!DEC$ENDIF
!
      END SUBROUTINE allreduce_r
!
!------------------------------------------------------------------------------
!
      SUBROUTINE allreducei_r(a,n)
!
!DEC$IF defined (mpi_flag)
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
!.....Input
      INTEGER :: n
!.....Input-output
      REAL(8) :: a(n)
!.....Local variables
      INTEGER :: i
      INTEGER :: ierr
      REAL(8) :: b(n)
!
      CALL MPI_ALLREDUCE(a,b,n,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!
      DO i=1,n
         a(i)=b(i)
      ENDDO
!
!DEC$ELSE
      IMPLICIT NONE
!
!.....Input
      INTEGER :: n
!.....Input-output
      REAL(8) :: a(n)
!DEC$ENDIF
!
      END SUBROUTINE allreducei_r
!
!------------------------------------------------------------------------------
!
      SUBROUTINE allreducei_r1(a)
!
!DEC$IF defined (mpi_flag)
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
!.....Input
      REAL(8) :: a
!.....Local variables
      INTEGER :: ierr
      REAL(8) :: b
!
      CALL MPI_ALLREDUCE(a,b,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      a=b
!
!DEC$ELSE
      IMPLICIT NONE
!
!.....Input
      REAL(8) :: a
!DEC$ENDIF
!
      END SUBROUTINE allreducei_r1
!
!------------------------------------------------------------------------------
!
      SUBROUTINE reduce_max_r1(a,b)
!
!DEC$IF defined (mpi_flag)
      USE Zcore        , ONLY: np
!
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!.....Input
      REAL(8) :: a
!.....Output
      REAL(8) :: b
!.....Local variables
      INTEGER root,ierr
!
      IF(np.gt.1) THEN
         root=0
         CALL MPI_REDUCE(a,b,1,MPI_DOUBLE_PRECISION,MPI_MAX, &
                         root,MPI_COMM_WORLD,ierr)
      ELSE
         b=a
      ENDIF
!DEC$ELSE
      IMPLICIT NONE
!
!.....Input
      REAL(8) :: a
!.....Output
      REAL(8) :: b
!
      b=a
!DEC$ENDIF
!
      END SUBROUTINE reduce_max_r1
!
!------------------------------------------------------------------------------
!
      SUBROUTINE reducei_max_r(a,n)
!
      USE Zcore  ,ONLY: myrank
!
!DEC$IF defined (mpi_flag)
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!.....Input
      INTEGER :: n
!.....Input-output
      REAL(8) :: a(n)
!.....Local variables
      INTEGER :: i,root,ierr
      REAL(8) :: b(n)
!
      root=0
      CALL MPI_REDUCE(a,b,n,MPI_DOUBLE_PRECISION,MPI_MAX, &
                      root,MPI_COMM_WORLD,ierr) 
!
      IF(myrank.eq.root) THEN
         DO i=1,n
            a(i)=b(i)
         ENDDO
      ENDIF
!
!DEC$ELSE
      IMPLICIT NONE
!
!.....Input
      INTEGER :: n
!.....Input-output
      REAL(8) :: a(n)
!DEC$ENDIF
!
      END SUBROUTINE reducei_max_r
!
!------------------------------------------------------------------------------
!
      SUBROUTINE allreducei_max_r(a,n)
!
!DEC$IF defined (mpi_flag)
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!.....Input
      INTEGER :: n
      REAL(8) :: a(n)
!.....Local variables
      INTEGER :: i
      INTEGER :: ierr
      REAL(8) ::  b(n)
!
      CALL MPI_ALLREDUCE(a,b,n,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
      DO i=1,n
         a(i)=b(i)
      ENDDO
!
!DEC$ELSE
      IMPLICIT NONE
!
!.....Input
      INTEGER :: n
      REAL(8) :: a(n)
!DEC$ENDIF
!
      END SUBROUTINE allreducei_max_r
!
!------------------------------------------------------------------------------
!
      SUBROUTINE allreducei_max_r1(a)
!
!DEC$IF defined (mpi_flag)
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!.....Input
      REAL(8) :: a
!.....Local variables
      INTEGER :: ierr
      REAL(8) :: b
!
      CALL MPI_ALLREDUCE(a,b,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
      a=b
!
!DEC$ELSE
      IMPLICIT NONE
!
!.....Input
      REAL(8) :: a
!DEC$ENDIF
!
      END SUBROUTINE allreducei_max_r1
!
!------------------------------------------------------------------------------
!
      SUBROUTINE allreduce_max_r1(a,b)
!
!DEC$IF defined (mpi_flag)
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!.....Input
      REAL(8) :: a
!.....Output
      REAL(8) :: b
!.....Local variables
      INTEGER :: ierr
!
      CALL MPI_ALLREDUCE(a,b,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
!
!DEC$ELSE
      IMPLICIT NONE
!
!.....Input
      REAL(8) :: a
!.....Output
      REAL(8) :: b
!
      b=a
!DEC$ENDIF
!
      END SUBROUTINE allreduce_max_r1
!
!------------------------------------------------------------------------------
!
      SUBROUTINE allreducei_min_r1(a)
!
!DEC$IF defined (mpi_flag)
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
!.....Input
      REAL(8) :: a
!.....Local variables
      INTEGER :: ierr
      REAL(8) :: b
!
      CALL MPI_ALLREDUCE(a,b,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
      a=b
!
!DEC$ELSE
      IMPLICIT NONE
!
!.....Input
      REAL(8) :: a
!DEC$ENDIF
!
      END SUBROUTINE allreducei_min_r1
!
!------------------------------------------------------------------------------
!
      SUBROUTINE allreduce_i(a,b,n)
!
!DEC$IF defined (mpi_flag)
      USE Zcore        , ONLY: np
!
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
!.....Input
      INTEGER :: n
      INTEGER :: a(n)
!.....Output
      INTEGER :: b(n)
!.....Local variables
      INTEGER :: i
      INTEGER :: ierr
!
      IF(np.gt.1) THEN
         CALL MPI_ALLREDUCE(a,b,n,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
      ELSE
         DO i=1,n
            b(i)=a(i)
         ENDDO
      ENDIF
!
!DEC$ELSE
      IMPLICIT NONE
!
!.....Input
      INTEGER :: n
      INTEGER :: a(n)
!.....Output
      INTEGER :: b(n)
!.....Local variables
      INTEGER :: i
!     
      DO i=1,n
         b(i)=a(i)
      ENDDO
!DEC$ENDIF
!
      END SUBROUTINE allreduce_i
!
!------------------------------------------------------------------------------
!
      SUBROUTINE allreduce_i1(a,b)
!
!DEC$IF defined (mpi_flag)
      USE Zcore        , ONLY: np
!
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
!.....Input
      INTEGER :: a
!.....Output
      INTEGER :: b
!.....Local variables
      INTEGER :: ierr
!
      IF(np.gt.1) THEN
         CALL MPI_ALLREDUCE(a,b,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
      ELSE
         b=a
      ENDIF
!
!DEC$ELSE
      IMPLICIT NONE
!
!.....Input
      INTEGER :: a
!.....Output
      INTEGER :: b
!     
      b=a
!DEC$ENDIF
!
      END SUBROUTINE allreduce_i1
!
!------------------------------------------------------------------------------
!
      SUBROUTINE reducei_i(a,n)
!
      USE Zcore           , ONLY: myrank
!
!DEC$IF defined (mpi_flag)
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
!.....Input
      INTEGER :: n
!.....Input-output
      INTEGER :: a(n)
!.....Local variables
      INTEGER :: i,root,ierr
      INTEGER :: b(n)
!
      root=0
      CALL MPI_REDUCE(a,b,n,MPI_INTEGER,MPI_SUM, &
                      root,MPI_COMM_WORLD,ierr)
!
      IF(myrank.eq.root) THEN
         DO i=1,n
            a(i)=b(i)
         ENDDO
      ENDIF
!
!DEC$ELSE
      IMPLICIT NONE
!
!.....Input
      INTEGER :: n
!.....Input-output
      INTEGER :: a(n)
!DEC$ENDIF
!
      END SUBROUTINE reducei_i
!
!------------------------------------------------------------------------------
!
      SUBROUTINE allreducei_i(a,n)
!
!DEC$IF defined (mpi_flag)
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
!.....Input
      INTEGER :: n
!.....Input-output
      INTEGER :: a(n)
!.....Local variables
      INTEGER :: i,ierr
      INTEGER :: b(n)
!
      CALL MPI_ALLREDUCE(a,b,n,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
!
      DO i=1,n
         a(i)=b(i)
      ENDDO
!
!DEC$ELSE
      IMPLICIT NONE
!
!.....Input
      INTEGER :: n
!.....Input-output
      INTEGER :: a(n)
!DEC$ENDIF
!
      END SUBROUTINE allreducei_i
!
!------------------------------------------------------------------------------
!
      SUBROUTINE allreducei_i1(a)
!
!DEC$IF defined (mpi_flag)
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
!.....Input-output
      INTEGER :: a
!.....Local variables
      INTEGER :: ierr
      INTEGER :: b
!
      CALL MPI_ALLREDUCE(a,b,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
      a=b
!
!DEC$ELSE
      IMPLICIT NONE
!
!.....Input-output
      INTEGER :: a
!DEC$ENDIF
!
      END SUBROUTINE allreducei_i1
!
!------------------------------------------------------------------------------
!
      SUBROUTINE allreduce_max_i(a,b,n)
!
!DEC$IF defined (mpi_flag)
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!.....Input
      INTEGER :: n
      INTEGER :: a(n)
!.....Output
      INTEGER :: b(n)
!.....Local variables
      INTEGER :: ierr
!
      CALL MPI_ALLREDUCE(a,b,n,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
!
!DEC$ELSE
      IMPLICIT NONE
!
!.....Input
      INTEGER :: n
      INTEGER :: a(n)
!.....Output
      INTEGER :: b(n)
!DEC$ENDIF
!
      END SUBROUTINE allreduce_max_i
!
!------------------------------------------------------------------------------
!
      SUBROUTINE allreducei_max_i1(a)
!
!DEC$IF defined (mpi_flag)
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!.....Input-output
      INTEGER :: a
!.....Local variables
      INTEGER :: ierr
      INTEGER :: b
!
      CALL MPI_ALLREDUCE(a,b,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
      a=b
!
!DEC$ELSE
      IMPLICIT NONE
!
!.....Input-output
      INTEGER :: a
!DEC$ENDIF
!
      END SUBROUTINE allreducei_max_i1
!
!------------------------------------------------------------------------------
!
      SUBROUTINE allreducei_l(a,n)
!
!DEC$IF defined (mpi_flag)
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
!.....Input
      INTEGER :: n
!.....Input-output
      LOGICAL :: a(n)
!.....Local variables
      INTEGER :: i,ierr
      LOGICAL :: b(n)
!
      CALL MPI_ALLREDUCE(a,b,n,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,ierr)
!
      DO i=1,n
         a(i)=b(i)
      ENDDO
!
!DEC$ELSE
      IMPLICIT NONE
!
!.....Input
      INTEGER :: n
!.....Input-output
      LOGICAL :: a(n)
!DEC$ENDIF
!
      END SUBROUTINE allreducei_l
!
!------------------------------------------------------------------------------
!
      SUBROUTINE allreducei_l1(a)
!
!DEC$IF defined (mpi_flag)
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
!.....Input-output
      LOGICAL :: a
!.....Local variables
      INTEGER :: ierr
      LOGICAL :: b
!
      CALL MPI_ALLREDUCE(a,b,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,ierr)
!
      a=b
!
!DEC$ELSE
      IMPLICIT NONE
!
!.....Input-output
      LOGICAL :: a
!DEC$ENDIF
!
      END SUBROUTINE allreducei_l1
!   
!------------------------------------------------------------------------------
!
      SUBROUTINE broadcast_i(a,n)
!
!DEC$IF defined (MCC)      
!DEC$ELSEIF defined (MCC_DLL)      
      !dec$ attributes dllexport :: bcast_i      
!DEC$ELSEIF defined (SPACE)          
      !dec$ attributes dllexport :: bcast_i
!DEC$ENDIF
!
!DEC$IF defined (mpi_flag)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
!
!.....Input
      INTEGER :: n
!.....Input-output
      INTEGER :: a(n)
!.....Local variables
      INTEGER :: ierr
!
      CALL MPI_BCAST(a,n,MPI_INTEGER,0,mpi_comm_world,ierr)
!DEC$ELSE
      IMPLICIT NONE
!
!.....Input
      INTEGER :: n
!.....Input-output
      INTEGER :: a(n)
!
!DEC$ENDIF
      END SUBROUTINE broadcast_i
!   
!------------------------------------------------------------------------------
!
      SUBROUTINE broadcast_i1(a)
!
!DEC$IF defined (MCC)      
!DEC$ELSEIF defined (MCC_DLL)      
      !dec$ attributes dllexport :: bcast_i1
!DEC$ELSEIF defined (SPACE)          
      !dec$ attributes dllexport :: bcast_i1
!DEC$ENDIF
!
!DEC$IF defined (mpi_flag)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
!
!.....Input-output
      INTEGER :: a
!.....Local variables
      INTEGER :: ierr
!
      CALL MPI_BCAST(a,1,MPI_INTEGER,0,mpi_comm_world,ierr)
!DEC$ELSE
      IMPLICIT NONE
!
!.....Input-output
      INTEGER :: a
!
!DEC$ENDIF
      END SUBROUTINE broadcast_i1
!   
!------------------------------------------------------------------------------
!
      SUBROUTINE broadcast_l1(a)
!
!DEC$IF defined (MCC)      
!DEC$ELSEIF defined (MCC_DLL)      
      !dec$ attributes dllexport :: bcast_i1
!DEC$ELSEIF defined (SPACE)          
      !dec$ attributes dllexport :: bcast_i1
!DEC$ENDIF
!
!DEC$IF defined (mpi_flag)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
!
!.....Input-output
      LOGICAL :: a
!.....Local variables
      INTEGER :: ierr
!
      CALL MPI_BCAST(a,1,MPI_LOGICAL,0,mpi_comm_world,ierr)
!DEC$ELSE
      IMPLICIT NONE
!
!.....Input-output
      LOGICAL :: a
!
!DEC$ENDIF
      END SUBROUTINE broadcast_l1
!   
!------------------------------------------------------------------------------
!
      SUBROUTINE broadcast_r(a,n)
!
!DEC$IF defined (MCC)      
!DEC$ELSEIF defined (MCC_DLL)      
      !dec$ attributes dllexport :: bcast_r      
!DEC$ELSEIF defined (SPACE)          
      !dec$ attributes dllexport :: bcast_r
!DEC$ENDIF
!
!DEC$IF defined (mpi_flag)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
!
!.....Input
      INTEGER :: n
!.....Input-output
      REAL(8) :: a(n)
!.....Local variables
      INTEGER :: ierr
!
      CALL mpi_bcast(a,n,MPI_DOUBLE_PRECISION,0,mpi_comm_world,ierr)
!DEC$ELSE
      IMPLICIT NONE
!
!.....Input
      INTEGER :: n
!.....Input-output
      REAL(8) :: a(n)
!
!DEC$ENDIF
      END SUBROUTINE broadcast_r
!   
!------------------------------------------------------------------------------
!
      SUBROUTINE broadcast_r1(a)
!
!DEC$IF defined (MCC)      
!DEC$ELSEIF defined (MCC_DLL)      
      !dec$ attributes dllexport :: bcast_r1
!DEC$ELSEIF defined (SPACE)          
      !dec$ attributes dllexport :: bcast_r1
!DEC$ENDIF
!
!DEC$IF defined (mpi_flag)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
!
!.....Input-output
      REAL(8) :: a
!.....Local variables
      INTEGER :: ierr
!
      CALL MPI_BCAST(a,1,MPI_DOUBLE_PRECISION,0,mpi_comm_world,ierr)
!
!DEC$ELSE
      IMPLICIT NONE
!
!.....Input-output
      REAL(8) :: a
!DEC$ENDIF
      END SUBROUTINE broadcast_r1
!
!------------------------------------------------------------------------------
!
      SUBROUTINE broadcast_rank_i1(a,rank)
!
!DEC$IF defined (mpi_flag)
!      
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
!.....Input
      INTEGER :: rank
!.....Input-output
      INTEGER :: a
!.....Local variables
      INTEGER :: ierr
!
      CALL MPI_BCAST(a,1,MPI_INTEGER,rank,mpi_comm_world,ierr)
!DEC$ELSE
!      
      IMPLICIT NONE
!
!.....Input
      INTEGER :: rank
!.....Input-output
      INTEGER :: a
!DEC$ENDIF
!
      END SUBROUTINE broadcast_rank_i1
!
!------------------------------------------------------------------------------
!
      SUBROUTINE broadcast_rank_r1(a,rank)
!
!DEC$IF defined (mpi_flag)
!      
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
!.....Input
      INTEGER :: rank
!.....Input-output
      REAL(8) :: a
!.....Local variables
      INTEGER :: ierr
!
      CALL MPI_BCAST(a,1,MPI_DOUBLE_PRECISION,rank,mpi_comm_world,ierr)
!DEC$ELSE
!      
      IMPLICIT NONE
!
!.....Input
      INTEGER :: rank
!.....Input-output
      REAL(8) :: a
!DEC$ENDIF
!
      END SUBROUTINE broadcast_rank_r1
!
!------------------------------------------------------------------------------
!
      SUBROUTINE finalize_mpi
!      
      IMPLICIT NONE
!
!DEC$IF defined (mpi_flag)
      INCLUDE 'mpif.h'
!
      INTEGER :: ierr
!
      CALL MPI_FINALIZE(ierr)
!DEC$ENDIF
!
      END SUBROUTINE finalize_mpi
!
!------------------------------------------------------------------------------
!
      SUBROUTINE barrier_mpi
!      
      IMPLICIT NONE
!DEC$IF defined (MCC)      
!DEC$ELSEIF defined (MCC_DLL)      
      !dec$ attributes dllexport :: barrier_mpi
!DEC$ELSEIF defined (SPACE)          
      !dec$ attributes dllexport :: barrier_mpi
!DEC$ENDIF
!      
!
!DEC$IF defined (mpi_flag)
      INCLUDE 'mpif.h'
!
      INTEGER ierr
!
      CALL mpi_barrier(mpi_comm_world,ierr)
!DEC$ENDIF
!
      END SUBROUTINE barrier_mpi
