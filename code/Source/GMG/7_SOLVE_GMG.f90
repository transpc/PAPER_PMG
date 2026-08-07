     SUBROUTINE SOLVE_GMG(icase)

! this routine solves pressure eq. by PMG
! icase = 1: for the new matrix (need to comput Galerkin condition RAP)
! icase = 2: for old matrix (doesnot need to comput RAP)
!

        USE MD_parameter, ONLY: mdf_matrix
        USE md_geometry, ONLY: nnode
        USE MD_matrix, ONLY: u, ut
        USE MD_MG_index, ONLY: isol_mg, isth   !,istiffness_MG
        USE MD_MPI, ONLY: myrank
        USE omp_lib

        IMPLICIT NONE

!DEC$IF defined (mpi_flag)
        INCLUDE 'mpif.h'
!DEC$ENDIF

        INTEGER(4) ierr, id, icase, i, j
        REAL(8) t1, t2

! - - - - - - - - - - - - - - - - - - - - - - - - !
!    SOLVE equation by PMG
!    Modified: April 2023.
!    Hybrid OpenMP+ MPI: Sept. 2023
! - - - - - - - - - - - - - - - - - - - - - - - - !

        ierr = 0
        id = 0

        !$omp PARALLEL DO
        DO i = 1, nnode
           ut(i) = u(i)
        END DO
        !$omp end PARALLEL DO

10      CONTINUE
        IF (ierr .EQ. 1) THEN
           mdf_matrix = 1 - mdf_matrix
        END IF

        IF (mdf_matrix .EQ. 1) THEN
           CALL Dig_mdf_matrix(icase)
        ELSE
           IF (ierr .EQ. 1) THEN
              CALL Dig_mdf_matrix_inv
           END IF
        END IF

! - - - - - - - - - - - - - - !
!     GALERKIN condition

!!DEC$IF defined (mpi_flag)
!      t1 = MPI_Wtime()
!!DEC$ENDIF

! - - - - - - - - - - - - - - - - -
        IF (icase == 1) THEN
           CALL stiffness_MG

! max_eig for Poly_smoothing

           IF (isth == 3) THEN

              CALL eig_value

           END IF
        END IF

! - - - - - - - - - - - - - - - -

!!DEC$IF defined (mpi_flag)
!      t2 = MPI_Wtime()

!         IF (myrank .eq. 0) WRITE (101, *) 'GAR out', t2-t1

!!DEC$ENDIF

! - - - - - - -
        IF (isol_mg .EQ. 1) THEN

           CALL SOLVER(ierr)

        ELSEIF (isol_mg == 2) THEN

           CALL SOLVER_NEW(ierr)

        ELSEIF (isol_mg .LE. 0) THEN

           CALL solve_pbcg_mg(ierr)
!            call solve_pbcg_ali(ierr)

        END IF
! - - - - - - -

        IF (ierr .EQ. 1) THEN

           !$omp PARALLEL DO
           DO i = 1, nnode
              u(i) = ut(i)
           END DO
           !$omp end PARALLEL DO

           id = id + 1
           IF (id .EQ. 2) THEN
              IF (myrank == 0) THEN
                 WRITE (999, *) 'divergence=> stop'
                 WRITE (*, *) 'divergence => stop'
              END IF

              STOP
           END IF

           GOTO 10
        END IF
!
        RETURN
     END

! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
     SUBROUTINE SOLVER(ierror)

! this routine for V-cycle PMG old version
!
        USE MD_MPI
        USE MD_MPI_MG
        USE MD_geometry, ONLY: nnode
        USE MD_matrix, ONLY: nnz, ia, ja, ju, au, u, b, ut
        USE MD_parameter, ONLY: crit, ndom
        USE MD_MG_matrix, ONLY: r, rt, rc, rs, e, et, es, &
                                iac, jac, juc, ias, jas, jus, &
                                iai, jai, iar, jar, &
                                auc, aus, Xrest, Xintp, nnzs
        USE MD_MG_index, ONLY: ncycle, iter_mg, nlevel, n_GC, &
                               relax, crit_1, maxit_1, isth
        USE MD_MG_coord, ONLY: ialv, nnods, ncolc, ncolf

!---------------------------------------------------------------------!
        IMPLICIT NONE
!DEC$IF defined (mpi_flag)
        INCLUDE 'mpif.h'
!DEC$ENDIF
        INTEGER(4) ierror
! --- temp
        INTEGER(4) icycle, i, np, j, k, i1, i2, Iter0
        INTEGER(4) ilv, ista, iend
        REAL(8) res, res0, res1
!
!DEC$IF defined (mpi_flag)
        INTEGER(4)::status(mpi_status_size), tag, ierr
!DEC$ENDIF

! ====================================================================!
! ---   Starting M-G Iteration method ---
! ====================================================================!
! --------
! S&R (u) * * * * * * * * * * * * * * * * * * * * * * * *

        np = ndom
        Iter0 = ITER_MG

        CALL send_receive(nnbd, nnode, spt, rpt, sintf, rintf, nbdom, u)
!
        CALL resi_normP(nnode, nintf, u, b, au, ja, ia, res1) ! res1 = (b-a*u)^2

!DEC$IF defined (mpi_flag)
!      CALL mpi_barrier(mpi_comm_world,ierr)
        res0 = 0.d0
        CALL mpi_allreduce(res1, res0, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
!DEC$ENDIF

        IF (np .GT. 1) THEN
           res0 = DSQRT(res0)
        ELSE
           res0 = DSQRT(res1)
        END IF

!     IF(myrank.eq.0) write(999,*)'res0=',res0

        rt = 0.d0
        rc = 0.d0
        r = 0.d0
        et = 0.d0

        DO 100 icycle = 1, ncycle

! ----
! 1. for finest level:
! ----PRE-SMOOTHING
           CALL smoothing_fine(Iter0, isth, ndom, relax, nintf, nnode, nnz, ia, ja, ju, au, &
                               u, b, nnbd, nbdom, spt, rpt, sintf, rintf)

!      call send_receive(nnbd,nnode,spt,rpt,sintf,rintf,nbdom,u)

! ---residual
           CALL resi_P(nnode, nintf, u, b, r, au, ja, ia)   ! r = b-a*u

! S&R (r)
           IF (nnbd .NE. 0) THEN
              CALL send_receive(nnbd, nnode, spt, rpt, sintf, rintf, nbdom, r)
           END IF
!
           rt(1:nnode) = r(1:nnode)

! 2: COARSER LEVELS:
           e = 0.d0
           DO ilv = 2, nlevel - 1
!     restriction: rc = R*rt;
!     Solve: A*e = rc
!     rt = rc-A*e

! restriction: rc = R*rt!

              ista = ialv(ilv) - nnode     ! for rc: from coarser grid
              iend = ista + iintf(ilv) - 1

              CALL matrix_vec_N(ista, iend, rc, rt, Xrest, iar, jar)

! smoothing: A*e=rc
!
              DO i = 1, Iter0
!
                 CALL Smooth_GS(1, ista, iend, rc, e, auc, iac, jac, juc)
!
                 CALL MD_S_R(ilv, ista, e)
!
              END DO

! rt = rc-A*e
              CALL residl(nnode, ista, iend, rt, rc, e, auc, iac, jac)
              ista = ista + nnode
              CALL MD_S_R(ilv, ista, rt)

           END DO

! 3: coarest level:
           ilv = nlevel
           ista = ialv(ilv) - nnode
           iend = ista + iintf(ilv) - 1     !ialv(ilv+1)-1 - nnode

           CALL matrix_vec_N(ista, iend, rc, rt, Xrest, iar, jar)

           rs(1:nintfs) = rc(ista:iend)
           es = 0.d0

           IF (n_GC .EQ. 0) THEN

              IF (iGS .EQ. 0) THEN
                 CALL PCG_Dig(maxit_1, nintfs, nnods, nnbds, nnzs, myrank, ias, jas, jus, &
                              spts, rpts, sintfs, rintfs, nbdoms, aus, rs, es, crit_1, ierror)
              ELSE

                 DO i = 1, nGS
                    CALL Relax_GSP(1, relax, nnods, nintfs, nnzs, ias, jas, jus, aus, es, rs)
                    CALL send_receive(nnbds, nnods, spts, rpts, sintfs, rintfs, nbdoms, es)
                 END DO
              END IF
!
           ELSE

              IF (icommu .EQ. 1) THEN
                 CALL SOLVE_GC(nintfs, nnods, maxit_1, crit_1, rs, es)
              ELSE
                 CALL SOLVE_GC_all(nintfs, nnods, maxit_1, crit_1, rs, es)
              END IF
!
           END IF

! 4: COARSER LEVELS:

           ilv = nlevel
           ista = ialv(ilv) - nnode
           iend = ialv(ilv + 1) - 1 - nnode

           e(ista:iend) = es(1:nnods)

           DO ilv = nlevel - 1, 2, -1

!     et = I*e
!     e = e + et
!     rc = A*e

              ista = ialv(ilv)                                  ! for et: from fine
              iend = ista + iintf(ilv) - 1     !ialv(ilv+1)-1

              CALL matrix_vec_N(ista, iend, et, e, Xintp, iai, jai)

              i1 = ista - nnode                                 ! for e.
              i2 = iend - nnode

              e(i1:i2) = e(i1:i2) + et(ista:iend)
              CALL MD_S_R(ilv, i1, e)
!
              !     Iter0 = Iter_mg
              !     IF(ilv.LE.3) Iter0 = Iter1

              DO i = 1, Iter0
                 CALL Smooth_GS(1, i1, i2, rc, e, auc, iac, jac, juc)
                 CALL MD_S_R(ilv, i1, e)
              END DO

           END DO

! 5: finest level

           ilv = 1
           ista = ialv(ilv)
           iend = ista + iintf(ilv) - 1 !ialv(ilv+1)-1

           CALL matrix_vec_N(ista, iend, et, e, Xintp, iai, jai)

           u(1:nintf) = u(1:nintf) + et(ista:iend)

! S&R (u)
           IF (nnbd .NE. 0) THEN
              CALL send_receive(nnbd, nnode, spt, rpt, sintf, rintf, nbdom, u)
           END IF

! POST SMOTHING

           CALL smoothing_fine(Iter0, isth, ndom, relax, nintf, nnode, nnz, ia, ja, ju, au, &
                               u, b, nnbd, nbdom, spt, rpt, sintf, rintf)

!      call send_receive(nnbd,nnode,spt,rpt,sintf,rintf,nbdom,u)

! ---residual-norm

           CALL resi_normP(nnode, nintf, u, b, au, ja, ia, res1)
!
! --
!DEC$IF defined (mpi_flag)
!      CALL mpi_barrier(mpi_comm_world,ierr)
           res = 0.d0
           CALL mpi_allreduce(res1, res, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
!DEC$ENDIF
           IF (np .GT. 1) THEN
              res = DSQRT(res)/(res0)
           ELSE
              res = DSQRT(res1)/(res0)
           END IF
!

           IF (res .lt. crit) GOTO 200

           IF (res .GT. 1.d6) THEN
              ierror = 1
              RETURN
           END IF
! ---
100     END DO   ! end main loop

        ierror = 1
        RETURN
! ---------------------
200     CONTINUE

        IF (myrank == 0) THEN
!      print*,'convergence, cycle=',icycle,res
           WRITE (16 + myrank, *) icycle
        END IF
! ---
        RETURN
     END

! = = = = = = = = = = = = = = = = = !
     SUBROUTINE matrix_vec(ista, iend, r, rt, Xr, ia, ja)
        IMPLICIT NONE
!
!-------------------------------------------------------------------------!
!                                                                         !
!  r = XR * rt                                                            !
!  input:                                                                 !
!     XR(n,nt)                                                            !
!     rt(nt) = array of length equal to the column dimension of matrix XR !
!     XR, ja, ia = input matrix in compressed sparse row format.          !
!  Output:                                                                !
!     r(n)  = real array of length n, containing the product r=XR*rt      !
! IT is noted that r and rt are long vector, each sub. only calculated    !
! r [ista,iend]                                                           !
!                                                                         !
!-------------------------------------------------------------------------!
!
        INTEGER(4) ista, iend
        INTEGER(4) ja(*)
        INTEGER(4) ia(*)
        REAL(8) Xr(*)
        REAL(8) r(*)
        REAL(8) rt(*)
! temp
        INTEGER i, k1, k2
!  ...
!  ...
        DO i = ista, iend
           k1 = ia(i)
           k2 = ia(i + 1) - 1
           r(i) = DOT_PRODUCT(Xr(k1:k2), rt(ja(k1:k2)))
        END DO
!=====
!=====
        RETURN
     END

! = = = = = = = = = = = = = = = = = !
     SUBROUTINE matrix_vec_N(ista, iend, y, x, a, ia, ja)

        USE omp_lib

        IMPLICIT NONE
!
!-------------------------------------------------------------------------!
!                                                                         !
!  r = XR * rt                                                            !
!  input:                                                                 !
!     XR(n,nt)                                                            !
!     rt(nt) = array of length equal to the column dimension of matrix XR !
!     XR, ja, ia = input matrix in compressed sparse row format.          !
!  Output:                                                                !
!     r(n)  = real array of length n, containing the product r=XR*rt      !
! IT is noted that r and rt are long vector, each sub. only calculated    !
! r [ista,iend]                                                           !
!                                                                         !
!-------------------------------------------------------------------------!
!
        INTEGER(4) ista, iend
        INTEGER(4) ja(*)
        INTEGER(4) ia(*)
        REAL(8) a(*)
        REAL(8) y(*)
        REAL(8) x(*)
! temp
        INTEGER i, k1, k2, k, l, m
        INTEGER j0, j1, j2, j3, j4, j5
        REAL(8) tmp
        REAL(8) t0, t1, t2, t3, t4, t5
!  ...
!  ...

        !$omp PARALLEL DO private(k1,k2,l,k,m,tmp,j0,j1,j2,j3,j4,j5,t0,t1,t2,t3,t4,t5)

        DO i = ista, iend
           k1 = ia(i)
           k2 = ia(i + 1)!-1

!   r(i) = DOT_PRODUCT( Xr(k1:k2),rt(ja(k1:k2)) )
!ENDDO

           l = k2 - k1 !ia(i+1)-ia(i)
           tmp = 0.d0
!
           IF (l .eq. 1) THEN
              k = k1! ia(i)
              j0 = ja(k)
              t0 = a(k)*x(j0)
              tmp = tmp + t0
           ELSEIF (l .eq. 2) THEN
              k = k1!ia(i)
              j0 = ja(k)
              j1 = ja(k + 1)
              t0 = a(k)*x(j0)
              t1 = a(k + 1)*x(j1)
              tmp = tmp + t0 + t1
           ELSEIF (l .eq. 3) THEN
              k = k1!ia(i)
              j0 = ja(k)
              j1 = ja(k + 1)
              j2 = ja(k + 2)
              t0 = a(k)*x(j0)
              t1 = a(k + 1)*x(j1)
              t2 = a(k + 2)*x(j2)
              tmp = tmp + t0 + t1 + t2
           ELSEIF (l .eq. 4) THEN
              k = k1! ia(i)
              j0 = ja(k)
              j1 = ja(k + 1)
              j2 = ja(k + 2)
              j3 = ja(k + 3)
              t0 = a(k)*x(j0)
              t1 = a(k + 1)*x(j1)
              t2 = a(k + 2)*x(j2)
              t3 = a(k + 3)*x(j3)
              tmp = tmp + t0 + t1 + t2 + t3
           ELSEIF (l .eq. 5) THEN
              k = k1! ia(i)
              j0 = ja(k)
              j1 = ja(k + 1)
              j2 = ja(k + 2)
              j3 = ja(k + 3)
              j4 = ja(k + 4)
              t0 = a(k)*x(j0)
              t1 = a(k + 1)*x(j1)
              t2 = a(k + 2)*x(j2)
              t3 = a(k + 3)*x(j3)
              t4 = a(k + 4)*x(j4)
              tmp = tmp + t0 + t1 + t2 + t3 + t4
           ELSEIF (l .eq. 6) THEN
              k = k1! ia(i)
              j0 = ja(k)
              j1 = ja(k + 1)
              j2 = ja(k + 2)
              j3 = ja(k + 3)
              j4 = ja(k + 4)
              j5 = ja(k + 5)
              t0 = a(k)*x(j0)
              t1 = a(k + 1)*x(j1)
              t2 = a(k + 2)*x(j2)
              t3 = a(k + 3)*x(j3)
              t4 = a(k + 4)*x(j4)
              t5 = a(k + 5)*x(j5)
              tmp = tmp + t0 + t1 + t2 + t3 + t4 + t5
           ELSE
              DO m = 0, l - 1
                 k = k1 + m ! ia(i  )+m
                 j0 = ja(k)
                 t0 = a(k)*x(j0)
                 tmp = tmp + t0
              END DO
           END IF
           y(i) = tmp
        END DO
        !$omp end PARALLEL DO
!=====
        RETURN
     END
! - - - - - - - - - -
! - - - - - - - - - - - - - - - - - - - - - - - - - - !
! = = = = = = = = = = = = = = = = = !
     SUBROUTINE matrix_vec_N_MPI(ista, iend, y, x, a, ia, ja)

!use omp_lib

        IMPLICIT NONE
!
!-------------------------------------------------------------------------!
!                                                                         !
!  r = XR * rt                                                            !
!  input:                                                                 !
!     XR(n,nt)                                                            !
!     rt(nt) = array of length equal to the column dimension of matrix XR !
!     XR, ja, ia = input matrix in compressed sparse row format.          !
!  Output:                                                                !
!     r(n)  = real array of length n, containing the product r=XR*rt      !
! IT is noted that r and rt are long vector, each sub. only calculated    !
! r [ista,iend]                                                           !
!                                                                         !
!-------------------------------------------------------------------------!
!
        INTEGER(4) ista, iend
        INTEGER(4) ja(*)
        INTEGER(4) ia(*)
        REAL(8) a(*)
        REAL(8) y(*)
        REAL(8) x(*)
! temp
        INTEGER i, k1, k2, k, l, m
        INTEGER j0, j1, j2, j3, j4, j5
        REAL(8) tmp
        REAL(8) t0, t1, t2, t3, t4, t5
!  ...
!  ...

!   !$omp PARALLEL DO private(k1,k2,l,k,m,tmp,j0,j1,j2,j3,j4,j5,t0,t1,t2,t3,t4,t5)

        DO i = ista, iend
           k1 = ia(i)
           k2 = ia(i + 1)!-1

!   r(i) = DOT_PRODUCT( Xr(k1:k2),rt(ja(k1:k2)) )
!ENDDO

           l = k2 - k1 !ia(i+1)-ia(i)
           tmp = 0.d0
!
           IF (l .eq. 1) THEN
              k = k1! ia(i)
              j0 = ja(k)
              t0 = a(k)*x(j0)
              tmp = tmp + t0
           ELSEIF (l .eq. 2) THEN
              k = k1!ia(i)
              j0 = ja(k)
              j1 = ja(k + 1)
              t0 = a(k)*x(j0)
              t1 = a(k + 1)*x(j1)
              tmp = tmp + t0 + t1
           ELSEIF (l .eq. 3) THEN
              k = k1!ia(i)
              j0 = ja(k)
              j1 = ja(k + 1)
              j2 = ja(k + 2)
              t0 = a(k)*x(j0)
              t1 = a(k + 1)*x(j1)
              t2 = a(k + 2)*x(j2)
              tmp = tmp + t0 + t1 + t2
           ELSEIF (l .eq. 4) THEN
              k = k1! ia(i)
              j0 = ja(k)
              j1 = ja(k + 1)
              j2 = ja(k + 2)
              j3 = ja(k + 3)
              t0 = a(k)*x(j0)
              t1 = a(k + 1)*x(j1)
              t2 = a(k + 2)*x(j2)
              t3 = a(k + 3)*x(j3)
              tmp = tmp + t0 + t1 + t2 + t3
           ELSEIF (l .eq. 5) THEN
              k = k1! ia(i)
              j0 = ja(k)
              j1 = ja(k + 1)
              j2 = ja(k + 2)
              j3 = ja(k + 3)
              j4 = ja(k + 4)
              t0 = a(k)*x(j0)
              t1 = a(k + 1)*x(j1)
              t2 = a(k + 2)*x(j2)
              t3 = a(k + 3)*x(j3)
              t4 = a(k + 4)*x(j4)
              tmp = tmp + t0 + t1 + t2 + t3 + t4
           ELSEIF (l .eq. 6) THEN
              k = k1! ia(i)
              j0 = ja(k)
              j1 = ja(k + 1)
              j2 = ja(k + 2)
              j3 = ja(k + 3)
              j4 = ja(k + 4)
              j5 = ja(k + 5)
              t0 = a(k)*x(j0)
              t1 = a(k + 1)*x(j1)
              t2 = a(k + 2)*x(j2)
              t3 = a(k + 3)*x(j3)
              t4 = a(k + 4)*x(j4)
              t5 = a(k + 5)*x(j5)
              tmp = tmp + t0 + t1 + t2 + t3 + t4 + t5
           ELSE
              DO m = 0, l - 1
                 k = k1 + m ! ia(i  )+m
                 j0 = ja(k)
                 t0 = a(k)*x(j0)
                 tmp = tmp + t0
              END DO
           END IF
           y(i) = tmp
        END DO
!      !$omp end PARALLEL DO
!=====
        RETURN
     END
! - - - - - - - - - -
! - - - - - - - - - - - - - - - - - - - - - - - - - - !

     SUBROUTINE residl(nnode, ista, iend, rt, b, u, Au, ia, ja)

        USE omp_lib

        IMPLICIT NONE
!
! rt = b-Au*u
! ---
        INTEGER(4) nnode, ista, iend
        INTEGER(4) ia(*)
        INTEGER(4) ja(*)
        REAL(8) b(*)
        REAL(8) u(*)
        REAL(8) au(*)
        REAL(8) rt(*)
! temp
        INTEGER(4) i, j1, j2, j, k
        REAL(8) temp
! ---

        !$omp PARALLEL DO private(j1,j2,temp,j)

        DO i = ista, iend
           j1 = ia(i)
           j2 = ia(i + 1) - 1
!
           temp = b(i)
           DO j = j1, j2
              temp = temp - au(j)*u(ja(j))
           END DO

! ---
!    temp = DOT_PRODUCT( au(j1:j2),u(ja(j1:j2)) )
           rt(i + nnode) = temp   ! b(i)-temp
        END DO

        !$omp end PARALLEL DO

        RETURN
     END

! = = = = = = = = = = = = = = = = = !
! - - - - - - - - - - - - - - - - - - - - - - - - - - !

     SUBROUTINE residl_MPI(nnode, ista, iend, rt, b, u, Au, ia, ja)

!use omp_lib

        IMPLICIT NONE
!
! rt = b-Au*u
! ---
        INTEGER(4) nnode, ista, iend
        INTEGER(4) ia(*)
        INTEGER(4) ja(*)
        REAL(8) b(*)
        REAL(8) u(*)
        REAL(8) au(*)
        REAL(8) rt(*)
! temp
        INTEGER(4) i, j1, j2, j, k
        REAL(8) temp
! ---

!      !$omp PARALLEL DO private(j1,j2,temp,j)

        DO i = ista, iend
           j1 = ia(i)
           j2 = ia(i + 1) - 1
!
           temp = b(i)
           DO j = j1, j2
              temp = temp - au(j)*u(ja(j))
           END DO

! ---
!    temp = DOT_PRODUCT( au(j1:j2),u(ja(j1:j2)) )
           rt(i + nnode) = temp   ! b(i)-temp
        END DO

!      !$omp end PARALLEL DO

        RETURN
     END

! = = = = = = = = = = = = = = = = = !

     SUBROUTINE Dig_mdf_matrix(icase)

! ---
        USE MD_matrix
        USE MD_MPI, ONLY: nintf
        USE MD_MG_matrix, ONLY: diagt
        USE omp_lib

!----------------------------------!
        IMPLICIT NONE

        INTEGER(4) icase
        INTEGER ie, i1, i2, nd, j, id, i

        REAL(8) xtmp

!
        IF (icase .EQ. 1) THEN

           !$omp PARALLEL DO private(i1,i2,j,xtmp,i)

           DO ie = 1, nintf
              i1 = ia(ie)
              i2 = ia(ie + 1) - 1
              j = ju(ie)
              xtmp = au(j)
!

              DO i = i1, i2
                 au(i) = au(i)/xtmp
              END DO
              b(ie) = b(ie)/xtmp
              diagr(ie) = 1.d0

           END DO

           !$omp end PARALLEL DO

        ELSE

           !$omp PARALLEL DO private(xtmp)

           DO ie = 1, nintf

              xtmp = diagt(ie)
              b(ie) = b(ie)/xtmp
           END DO

           !$omp end PARALLEL DO

        END IF

        !
        RETURN
     END

! = = = = = = = = = = = = = = = = = !

     SUBROUTINE Dig_mdf_matrix_inv

! ---
        USE MD_matrix
        USE MD_MPI, ONLY: nintf
        USE MD_MG_matrix, ONLY: diagt
        USE omp_lib

!----------------------------------!
        IMPLICIT NONE

        INTEGER ie, i1, i2, nd, j, id, i

        REAL(8) xtmp
!
        !$omp PARALLEL DO private(i1,i2,xtmp,i)

        DO ie = 1, nintf
           i1 = ia(ie)
           i2 = ia(ie + 1) - 1
           !
           xtmp = diagt(ie)

           DO i = i1, i2
              au(i) = au(i)*xtmp
           END DO

           b(ie) = b(ie)*xtmp
           diagr(ie) = 1.d0/xtmp

        END DO

        !$omp end PARALLEL DO

        !

        RETURN
     END

! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
     SUBROUTINE SOLVER_NEW(ierror)
!
        USE MD_MPI, ONLY: nintf, myrank
        USE MD_MPI_MG, ONLY: icommu, iGS, nGS, nintfs, nnbds, spts, rpts, sintfs, &
                             rintfs, nbdoms, iintf
        USE MD_geometry, ONLY: nnode
        USE MD_matrix, ONLY: nnz, ia, ja, ju, au, u, b, alu
        USE MD_parameter, ONLY: crit, ndom
        USE MD_MG_matrix, ONLY: r, rt, rc, rs, e, et, es, &
                                iac, jac, juc, ias, jas, jus, &
                                iai, jai, iar, jar, &
                                auc, aus, Xrest, Xintp, nnzs, diagrc
        USE MD_MG_index, ONLY: ncycle, iter_mg, nlevel, n_GC, &
                               relax, crit_1, maxit_1, isth, id_GS_sym, itergs
        USE MD_MG_coord, ONLY: ialv, nnods, ncolc, ncolf

        USE MD_MPI_ARP, ONLY: nnbdA, sptA, rptA, sintfA, rintfA, nbdomA, &
                              nnbdR, sptR, rptR, sintfR, rintfR, nbdomR

        USE omp_lib
!---------------------------------------------------------------------!
        IMPLICIT NONE
!DEC$IF defined (mpi_flag)
        INCLUDE 'mpif.h'
!DEC$ENDIF
        INTEGER(4) ierror
! --- temp
        INTEGER(4) icycle, i, np, j, k, i1, i2, Iter0, iter1, id
        INTEGER(4) ilv, ista, iend
        REAL(8) res, res0, res1, t1, t2, t3
!
!DEC$IF defined (mpi_flag)
        INTEGER(4)::status(mpi_status_size), tag, ierr
!DEC$ENDIF

! ====================================================================!
! ---   Starting M-G Iteration method ---
! ====================================================================!
! --------
        np = ndom
!      Iter0 = ITER_MG
!      iter1=1
        Iter0 = itergs(1)

! for ILU smoothing:
        IF (isth == 2) THEN
           alu = au
           CALL ilupcp(nintf, nnode, nnz, ia, ja, ju, alu)
! test:
!        alu = 0.d0
!      CALL ilupcp_new(nintf,nnode,nnz,ia,ja,ju,au,alu)

        END IF

! NEW for A - - - - - - -
        IF (nnbdA .NE. 0) THEN
           CALL send_receive(nnbdA, nnode, sptA, rptA, sintfA, rintfA, nbdomA, u)
        END IF
!
        CALL resi_normP(nnode, nintf, u, b, au, ja, ia, res1) ! res1 = (b-a*u)^2

!DEC$IF defined (mpi_flag)
!      CALL mpi_barrier(mpi_comm_world,ierr)
        res0 = 0.d0
        CALL mpi_allreduce(res1, res0, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
!DEC$ENDIF

        IF (np .GT. 1) THEN
           res0 = DSQRT(res0)
        ELSE
           res0 = DSQRT(res1)
        END IF

!     IF(myrank.eq.0) write(999,*)'res0=',res0

        IF (res0 .LE. 1.d-16) THEN
           IF (myrank .EQ. 0) WRITE (999, *) 'ncycle=', 0
           RETURN
        END IF

!      rt = 0.d0
!      rc = 0.d0
!      r = 0.d0
!      et = 0.d0
!
        DO 100 icycle = 1, ncycle

!
           Iter0 = itergs(1)
! ----
! 1. for finest level:
! ----PRE-SMOOTHING
! NEW for A
           CALL smoothing_fine(Iter0, isth, ndom, relax, nintf, nnode, nnz, ia, ja, ju, au, &
                               u, b, nnbdA, nbdomA, sptA, rptA, sintfA, rintfA)

!      call send_receive(nnbd,nnode,spt,rpt,sintf,rintf,nbdom,u)

! ---residual
           CALL resi_P(nnode, nintf, u, b, r, au, ja, ia)   ! r = b-a*u

! S&R (r)

! NEW for R
           IF (nnbdR .NE. 0) THEN
              CALL send_receive(nnbdR, nnode, sptR, rptR, sintfR, rintfR, nbdomR, r)
           END IF
!
           !$omp PARALLEL DO
           DO i = 1, nnode
              rt(i) = r(i)
           END DO
           !$omp end PARALLEL DO

! 2: COARSER LEVELS: = = = = = = = = = = = = = = = = = = = = =

           !$omp PARALLEL DO
           DO i = 1, ncolc
              e(i) = 0.d0
           END DO
           !$omp end PARALLEL DO

           DO ilv = 2, nlevel - 1
!     restriction: rc = R*rt;
!     Solve: A*e = rc
!     rt = rc-A*e

              Iter0 = itergs(ilv)

! restriction: rc = R*rt!

              ista = ialv(ilv) - nnode     ! for rc: from coarser grid
              iend = ista + iintf(ilv) - 1

              CALL matrix_vec_N(ista, iend, rc, rt, Xrest, iar, jar)

! smoothing: A*e=rc
!

!      IF(ilv.EQ.2) THEN
!          Iter0 = Iter1
!      ELSE
!          Iter0 = ITER_MG
!      ENDIF
!

              DO i = 1, Iter0
!
                 IF ((MOD(i, 2) .EQ. 1) .OR. (id_GS_sym .EQ. 0)) THEN

                    CALL Smooth_GS2(1, ista, iend, rc, e, auc, iac, jac, juc, diagrc)

                 ELSE
                    CALL Smooth_GS_BW(1, ista, iend, rc, e, auc, iac, jac, juc)
                 END IF

!
! NEW for A
                 id = 1
                 CALL MD_S_R_NEW(id, ilv, ista, e)
!
              END DO

! rt = rc-A*e
              CALL residl(nnode, ista, iend, rt, rc, e, auc, iac, jac)
              ista = ista + nnode
! NEW for R
              id = 2
              CALL MD_S_R_NEW(id, ilv, ista, rt)

           END DO

! 3: coarest level:
           ilv = nlevel
           ista = ialv(ilv) - nnode
           iend = ista + iintf(ilv) - 1     !ialv(ilv+1)-1 - nnode

           CALL matrix_vec_N(ista, iend, rc, rt, Xrest, iar, jar)

           j = ista - 1
           !$omp PARALLEL DO
           DO i = 1, nintfs
              rs(i) = rc(i + j)
              es(i) = 0.d0
           END DO
           !$omp END PARALLEL DO

           IF (n_GC .EQ. 0) THEN

              IF (iGS .EQ. 0) THEN
                 CALL PCG_Dig(maxit_1, nintfs, nnods, nnbds, nnzs, myrank, ias, jas, jus, &
                              spts, rpts, sintfs, rintfs, nbdoms, aus, rs, es, crit_1, ierror)
              ELSE

                 DO i = 1, nGS
                    CALL Relax_GSP(1, relax, nnods, nintfs, nnzs, ias, jas, jus, aus, es, rs)
                    CALL send_receive(nnbds, nnods, spts, rpts, sintfs, rintfs, nbdoms, es)
                 END DO
              END IF
!
           ELSE

              IF (icommu .EQ. 1) THEN
                 CALL SOLVE_GC(nintfs, nnods, maxit_1, crit_1, rs, es)
              ELSE
                 CALL SOLVE_GC_all(nintfs, nnods, maxit_1, crit_1, rs, es)
              END IF
!
           END IF

! 4: COARSER LEVELS:

           ilv = nlevel
           ista = ialv(ilv) - nnode
           iend = ialv(ilv + 1) - 1 - nnode

           j = ista - 1
           !$omp PARALLEL DO
           DO i = 1, nnods
              e(i + j) = es(i)
           END DO
           !$omp  END PARALLEL DO

           DO ilv = nlevel - 1, 2, -1

!     et = I*e
!     e = e + et
!     rc = A*e

              Iter0 = itergs(ilv)
!
              ista = ialv(ilv)                                  ! for et: from fine
              iend = ista + iintf(ilv) - 1     !ialv(ilv+1)-1

              CALL matrix_vec_N(ista, iend, et, e, Xintp, iai, jai)

              i1 = ista - nnode                                 ! for e.
              i2 = iend - nnode

              !$omp PARALLEL DO
              DO i = i1, i2
                 e(i) = e(i) + et(i + nnode)
              END DO
              !$omp END PARALLEL DO

!
!      IF(ilv.EQ.2) THEN
!          Iter0 = Iter1
!      ELSE
!          Iter0 = ITER_MG
!      ENDIF
!
! NEW for A
              id = 1

              DO i = 1, Iter0
                 CALL MD_S_R_NEW(id, ilv, i1, e)

                 IF ((MOD(i, 2) .EQ. 1) .OR. (id_GS_sym .EQ. 0)) THEN
                    CALL Smooth_GS2(1, i1, i2, rc, e, auc, iac, jac, juc, diagrc)
                 ELSE
                    CALL Smooth_GS_BW(1, i1, i2, rc, e, auc, iac, jac, juc)
                 END IF

              END DO

! NEW for P

              id = 3

              CALL MD_S_R_NEW(id, ilv, i1, e)

           END DO

! 5: finest level

           ilv = 1
           ista = ialv(ilv)
           iend = ista + iintf(ilv) - 1 !ialv(ilv+1)-1

           CALL matrix_vec_N(ista, iend, et, e, Xintp, iai, jai)

           j = ista - 1
           !$omp PARALLEL DO
           DO i = 1, nintf
              u(i) = u(i) + et(i + j)
           END DO
           !$omp END PARALLEL DO

! S&R (u)
! NEW for A
           IF (nnbdA .NE. 0) THEN
              CALL send_receive(nnbdA, nnode, sptA, rptA, sintfA, rintfA, nbdomA, u)
           END IF

! POST SMOTHING
           Iter0 = itergs(1)
! NEW For A
           CALL smoothing_fine(Iter0, isth, ndom, relax, nintf, nnode, nnz, ia, ja, ju, au, &
                               u, b, nnbdA, nbdomA, sptA, rptA, sintfA, rintfA)

!      call send_receive(nnbd,nnode,spt,rpt,sintf,rintf,nbdom,u)

! ---residual-norm

           CALL resi_normP(nnode, nintf, u, b, au, ja, ia, res1)
!
! --
!DEC$IF defined (mpi_flag)
!      CALL mpi_barrier(mpi_comm_world,ierr)
           res = 0.d0
           CALL mpi_allreduce(res1, res, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
!DEC$ENDIF
           IF (np .GT. 1) THEN
              res = DSQRT(res)/(res0)
           ELSE
              res = DSQRT(res1)/(res0)
           END IF
!
!      if(myrank==0) then
!      print*,'cycle=',icycle,res
!      ENDIF

           IF (res .lt. crit) GOTO 200

           IF (res .GT. 1.d6) THEN        ! divergence
              ierror = 1
              RETURN
           END IF
!
! ---
100     END DO   ! end main loop

        ierror = 1
        RETURN
! ---------------------
200     CONTINUE

        IF (myrank == 0) THEN
!      print*,'convergence, cycle=',icycle,res
           WRITE (16 + myrank, *) icycle
        END IF
! ---
        RETURN
     END

! = = = = = = = = = = = = = = = = = !

! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
     SUBROUTINE SOLVER_NEW_MPI(ierror)
!
        USE MD_MPI, ONLY: nintf, myrank
        USE MD_MPI_MG, ONLY: icommu, iGS, nGS, nintfs, nnbds, spts, rpts, sintfs, &
                             rintfs, nbdoms, iintf
        USE MD_geometry, ONLY: nnode
        USE MD_matrix, ONLY: nnz, ia, ja, ju, au, u, b, alu
        USE MD_parameter, ONLY: crit, ndom
        USE MD_MG_matrix, ONLY: r, rt, rc, rs, e, et, es, &
                                iac, jac, juc, ias, jas, jus, &
                                iai, jai, iar, jar, &
                                auc, aus, Xrest, Xintp, nnzs, diagrc
        USE MD_MG_index, ONLY: ncycle, iter_mg, nlevel, n_GC, &
                               relax, crit_1, maxit_1, isth, id_GS_sym, itergs
        USE MD_MG_coord, ONLY: ialv, nnods, ncolc, ncolf

        USE MD_MPI_ARP, ONLY: nnbdA, sptA, rptA, sintfA, rintfA, nbdomA, &
                              nnbdR, sptR, rptR, sintfR, rintfR, nbdomR

!     use omp_lib
!---------------------------------------------------------------------!
        IMPLICIT NONE
!DEC$IF defined (mpi_flag)
        INCLUDE 'mpif.h'
!DEC$ENDIF
        INTEGER(4) ierror
! --- temp
        INTEGER(4) icycle, i, np, j, k, i1, i2, Iter0, iter1, id
        INTEGER(4) ilv, ista, iend
        REAL(8) res, res0, res1, t1, t2, t3
!
!DEC$IF defined (mpi_flag)
        INTEGER(4)::status(mpi_status_size), tag, ierr
!DEC$ENDIF

! ====================================================================!
! ---   Starting M-G Iteration method ---
! ====================================================================!
! --------
        np = ndom
!      Iter0 = ITER_MG
!      iter1=1
        Iter0 = itergs(1)

! for ILU smoothing:
        IF (isth == 2) THEN
           alu = au
           CALL ilupcp(nintf, nnode, nnz, ia, ja, ju, alu)
        END IF

! NEW for A - - - - - - -
        IF (nnbdA .NE. 0) THEN
           CALL send_receive(nnbdA, nnode, sptA, rptA, sintfA, rintfA, nbdomA, u)
        END IF
!
        CALL resi_normP_MPI(nnode, nintf, u, b, au, ja, ia, res1) ! res1 = (b-a*u)^2

!DEC$IF defined (mpi_flag)
!      CALL mpi_barrier(mpi_comm_world,ierr)
        res0 = 0.d0
        CALL mpi_allreduce(res1, res0, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
!DEC$ENDIF

        IF (np .GT. 1) THEN
           res0 = DSQRT(res0)
        ELSE
           res0 = DSQRT(res1)
        END IF

!     IF(myrank.eq.0) write(999,*)'res0=',res0

        IF (res0 .LE. 1.d-20) THEN
           IF (myrank .EQ. 0) WRITE (999, *) 'ncycle=', 0
           RETURN
        END IF

!      rt = 0.d0
!      rc = 0.d0
!      r = 0.d0
!      et = 0.d0
!
        DO 100 icycle = 1, ncycle

!
           Iter0 = itergs(1)
! ----
! 1. for finest level:
! ----PRE-SMOOTHING
! NEW for A
           CALL smoothing_fine(Iter0, isth, ndom, relax, nintf, nnode, nnz, ia, ja, ju, au, &
                               u, b, nnbdA, nbdomA, sptA, rptA, sintfA, rintfA)

!      call send_receive(nnbd,nnode,spt,rpt,sintf,rintf,nbdom,u)

! ---residual
           CALL resi_P_MPI(nnode, nintf, u, b, r, au, ja, ia)   ! r = b-a*u

! S&R (r)

! NEW for R
           IF (nnbdR .NE. 0) THEN
              CALL send_receive(nnbdR, nnode, sptR, rptR, sintfR, rintfR, nbdomR, r)
           END IF
!
!    !$omp PARALLEL DO
           DO i = 1, nnode
              rt(i) = r(i)
           END DO
!    !$omp end PARALLEL DO

! 2: COARSER LEVELS: = = = = = = = = = = = = = = = = = = = = =

!    !$omp PARALLEL DO
           DO i = 1, ncolc
              e(i) = 0.d0
           END DO
!    !$omp end PARALLEL DO

           DO ilv = 2, nlevel - 1
!     restriction: rc = R*rt;
!     Solve: A*e = rc
!     rt = rc-A*e

              Iter0 = itergs(ilv)

! restriction: rc = R*rt!

              ista = ialv(ilv) - nnode     ! for rc: from coarser grid
              iend = ista + iintf(ilv) - 1

              CALL matrix_vec_N_MPI(ista, iend, rc, rt, Xrest, iar, jar)

! smoothing: A*e=rc
!

!      IF(ilv.EQ.2) THEN
!          Iter0 = Iter1
!      ELSE
!          Iter0 = ITER_MG
!      ENDIF
!

              DO i = 1, Iter0
!
                 IF ((MOD(i, 2) .EQ. 1) .OR. (id_GS_sym .EQ. 0)) THEN

                    CALL Smooth_GS2_MPI(1, ista, iend, rc, e, auc, iac, jac, juc, diagrc)

                 ELSE
                    CALL Smooth_GS_BW(1, ista, iend, rc, e, auc, iac, jac, juc)
                 END IF

!
! NEW for A
                 id = 1
                 CALL MD_S_R_NEW(id, ilv, ista, e)
!
              END DO

! rt = rc-A*e
              CALL residl_MPI(nnode, ista, iend, rt, rc, e, auc, iac, jac)
              ista = ista + nnode
! NEW for R
              id = 2
              CALL MD_S_R_NEW(id, ilv, ista, rt)

           END DO

! 3: coarest level:
           ilv = nlevel
           ista = ialv(ilv) - nnode
           iend = ista + iintf(ilv) - 1     !ialv(ilv+1)-1 - nnode

           CALL matrix_vec_N_MPI(ista, iend, rc, rt, Xrest, iar, jar)

           j = ista - 1
!    !$omp PARALLEL DO
           DO i = 1, nintfs
              rs(i) = rc(i + j)
              es(i) = 0.d0
           END DO
!    !$omp END PARALLEL DO

           IF (n_GC .EQ. 0) THEN

              IF (iGS .EQ. 0) THEN
                 CALL PCG_Dig(maxit_1, nintfs, nnods, nnbds, nnzs, myrank, ias, jas, jus, &
                              spts, rpts, sintfs, rintfs, nbdoms, aus, rs, es, crit_1, ierror)
              ELSE

                 DO i = 1, nGS
                    CALL Relax_GSP(1, relax, nnods, nintfs, nnzs, ias, jas, jus, aus, es, rs)
                    CALL send_receive(nnbds, nnods, spts, rpts, sintfs, rintfs, nbdoms, es)
                 END DO
              END IF
!
           ELSE

              IF (icommu .EQ. 1) THEN
                 CALL SOLVE_GC(nintfs, nnods, maxit_1, crit_1, rs, es)
              ELSE
                 CALL SOLVE_GC_all_MPI(nintfs, nnods, maxit_1, crit_1, rs, es)
              END IF
!
           END IF

! 4: COARSER LEVELS:

           ilv = nlevel
           ista = ialv(ilv) - nnode
           iend = ialv(ilv + 1) - 1 - nnode

           j = ista - 1
!    !$omp PARALLEL DO
           DO i = 1, nnods
              e(i + j) = es(i)
           END DO
!    !$omp  END PARALLEL DO

           DO ilv = nlevel - 1, 2, -1

!     et = I*e
!     e = e + et
!     rc = A*e

              Iter0 = itergs(ilv)
!
              ista = ialv(ilv)                                  ! for et: from fine
              iend = ista + iintf(ilv) - 1     !ialv(ilv+1)-1

              CALL matrix_vec_N_MPI(ista, iend, et, e, Xintp, iai, jai)

              i1 = ista - nnode                                 ! for e.
              i2 = iend - nnode

!    !$omp PARALLEL DO
              DO i = i1, i2
                 e(i) = e(i) + et(i + nnode)
              END DO
!    !$omp END PARALLEL DO

!
!      IF(ilv.EQ.2) THEN
!          Iter0 = Iter1
!      ELSE
!          Iter0 = ITER_MG
!      ENDIF
!
! NEW for A
              id = 1

              DO i = 1, Iter0
                 CALL MD_S_R_NEW(id, ilv, i1, e)

                 IF ((MOD(i, 2) .EQ. 1) .OR. (id_GS_sym .EQ. 0)) THEN
                    CALL Smooth_GS2_MPI(1, i1, i2, rc, e, auc, iac, jac, juc, diagrc)
                 ELSE
                    CALL Smooth_GS_BW(1, i1, i2, rc, e, auc, iac, jac, juc)
                 END IF

              END DO

! NEW for P

              id = 3

              CALL MD_S_R_NEW(id, ilv, i1, e)

           END DO

! 5: finest level

           ilv = 1
           ista = ialv(ilv)
           iend = ista + iintf(ilv) - 1 !ialv(ilv+1)-1

           CALL matrix_vec_N_MPI(ista, iend, et, e, Xintp, iai, jai)

           j = ista - 1
!    !$omp PARALLEL DO
           DO i = 1, nintf
              u(i) = u(i) + et(i + j)
           END DO
           !   !$omp END PARALLEL DO

! S&R (u)
! NEW for A
           IF (nnbdA .NE. 0) THEN
              CALL send_receive(nnbdA, nnode, sptA, rptA, sintfA, rintfA, nbdomA, u)
           END IF

! POST SMOTHING
           Iter0 = itergs(1)
! NEW For A
           CALL smoothing_fine(Iter0, isth, ndom, relax, nintf, nnode, nnz, ia, ja, ju, au, &
                               u, b, nnbdA, nbdomA, sptA, rptA, sintfA, rintfA)

!      call send_receive(nnbd,nnode,spt,rpt,sintf,rintf,nbdom,u)

! ---residual-norm

           CALL resi_normP_MPI(nnode, nintf, u, b, au, ja, ia, res1)
!
! --
!DEC$IF defined (mpi_flag)
!      CALL mpi_barrier(mpi_comm_world,ierr)
           res = 0.d0
           CALL mpi_allreduce(res1, res, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
!DEC$ENDIF
           IF (np .GT. 1) THEN
              res = DSQRT(res)/(res0)
           ELSE
              res = DSQRT(res1)/(res0)
           END IF
!
!      if(myrank==0) then
!      print*,'cycle=',icycle,res
!      ENDIF

           IF (res .lt. crit) GOTO 200

           IF (res .GT. 1.d6) THEN        ! divergence
              ierror = 1
              RETURN
           END IF
!
! ---
100     END DO   ! end main loop

        ierror = 1
        RETURN
! ---------------------
200     CONTINUE

        IF (myrank == 0) THEN
!      print*,'convergence, cycle=',icycle,res
           WRITE (16 + myrank, *) icycle
        END IF
! ---
        RETURN
     END

! = = = = = = = = = = = = = = = = = !

