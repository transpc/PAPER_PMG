        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun  3 00:45:09 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PARDISO_SOLVE__genmod
          INTERFACE 
            SUBROUTINE PARDISO_SOLVE(N,IA,JA,A,X,B)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: IA(1)
              INTEGER(KIND=4) :: JA(1)
              REAL(KIND=8) :: A(1)
              REAL(KIND=8) :: X(1)
              REAL(KIND=8) :: B(1)
            END SUBROUTINE PARDISO_SOLVE
          END INTERFACE 
        END MODULE PARDISO_SOLVE__genmod
