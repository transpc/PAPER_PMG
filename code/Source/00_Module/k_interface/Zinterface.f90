      MODULE Zinterface
         INTERFACE
            SUBROUTINE communicate_1d(a1,a2,a3,a4,a5,a6,a7,a8)
              USE Zmpi     , ONLY: ncell_fp
              REAL(8),INTENT(INOUT) :: a1(ncell_fp)
              REAL(8),DIMENSION(ncell_fp),OPTIONAL :: a2,a3,a4,a5,a6,a7,a8
            END SUBROUTINE communicate_1d
         END INTERFACE
         INTERFACE
            SUBROUTINE communicate_2d(a1,a2,a3)
              USE Zmpi     , ONLY: ncell_fp
              USE Zparam   , ONLY: ndim
              REAL(8),INTENT(INOUT) :: a1(ncell_fp,ndim)
              REAL(8),DIMENSION(ncell_fp,ndim),OPTIONAL :: a2,a3
            END SUBROUTINE communicate_2d
         END INTERFACE
         INTERFACE
            SUBROUTINE communicate_3d(a1,a2,a3)
              USE Zmpi     , ONLY: ncell_fp
              USE Zparam   , ONLY: ndim
              REAL(8),INTENT(INOUT) :: a1(ncell_fp,ndim*ndim)
              REAL(8),DIMENSION(ncell_fp,ndim*ndim),OPTIONAL :: a2,a3
            END SUBROUTINE communicate_3d
         END INTERFACE
         INTERFACE
            SUBROUTINE communicate_1d_int(a1,a2,a3)
              USE Zmpi     , ONLY: ncell_fp
              INTEGER,INTENT(INOUT) :: a1(ncell_fp)
              INTEGER,DIMENSION(ncell_fp),OPTIONAL :: a2,a3
            END SUBROUTINE communicate_1d_int
         END INTERFACE
         INTERFACE
            SUBROUTINE communicate_rv_2d(a1,a2)
              USE Zrv_mpi      , ONLY: ncell_fuel_rod_p
              USE Zrv_hts_2d   , ONLY: nr_2d
              REAL(8),INTENT(INOUT) :: a1(ncell_fuel_rod_p,nr_2d)
              REAL(8),DIMENSION(ncell_fuel_rod_p,nr_2d),OPTIONAL :: a2
            END SUBROUTINE communicate_rv_2d
         END INTERFACE
         INTERFACE
            SUBROUTINE sum_nf(zero,sym, &
                              s1i,s1,   &
                              s2i,s2,   &
                              s3i,s3,   &
                              s4i,s4,   &
                              s5i,s5,   &
                              s6i,s6,   &
                              s7i,s7,   &
                              s8i,s8)
            USE Zzone        , ONLY: ncell_fluid
            USE Znum_cell    , ONLY: lens
            INTEGER,INTENT(IN) :: zero,sym
            REAL(8),INTENT(IN) :: s1i(lens)
            REAL(8),INTENT(OUT) :: s1(ncell_fluid)
            REAL(8),DIMENSION(lens),OPTIONAL :: s2i,s3i,s4i,s5i,s6i,s7i,s8i
            REAL(8),DIMENSION(ncell_fluid),OPTIONAL :: s2,s3,s4,s5,s6,s7,s8
            END SUBROUTINE sum_nf
         END INTERFACE
         INTERFACE
            SUBROUTINE sum_nf_ndim(zero,sym,ncell, &
                                   s1i,s1,       &
                                   s2i,s2,       &
                                   s3i,s3,       &
                                   s4i,s4)
            USE Zzone        , ONLY: ncell_fluid
            USE Zparam       , ONLY: ndim
            USE Znum_cell    , ONLY: lens
            INTEGER,INTENT(IN) :: zero,sym
            INTEGER :: ncell
            REAL(8),INTENT(IN) :: s1i(lens,ndim)
            REAL(8),DIMENSION(lens,ndim),OPTIONAL :: s2i,s3i,s4i
            REAL(8),INTENT(OUT) :: s1(ncell,ndim)
            REAL(8),DIMENSION(ncell,ndim),OPTIONAL :: s2,s3,s4
            END SUBROUTINE sum_nf_ndim
         END INTERFACE
         INTERFACE
            SUBROUTINE pbcg_ilu(eps,maxiter,ncell,ncell_pad,maxmt_pad,maxmt_lu0,maxmt_lu1, &
                                diag,                                                      &
                                diag_lu,alu0,alu1,ia0,ia1,ja0,ja1,                         &
                                iap,jap,ap,jaar,iaa,ngroup,nbgroup,                        &
                                lev_typet,perm_r,                                          &
                                neq,arhsu,solu,izone,isPSolve)
            INTEGER :: maxiter,ncell,ncell_pad,maxmt_pad,maxmt_lu0,maxmt_lu1
            INTEGER :: neq,lev_typet
            INTEGER :: izone
            INTEGER :: ia0(ncell+1),ia1(ncell+1),ja0(maxmt_lu0),ja1(maxmt_lu1)
            INTEGER :: iap(2,ngroup+1),jap(maxmt_pad),jaar(ncell),iaa(2,ngroup+1),ngroup,nbgroup(3,ngroup)
            INTEGER :: perm_r(ncell)
            REAL(8) :: diag(ncell)
            REAL*8  :: diag_lu(ncell),alu0(maxmt_lu0),alu1(maxmt_lu1)
            REAL(8) :: ap(maxmt_pad)
            REAL(8) :: arhsu(ncell)
            REAL(8) :: solu(neq)
            LOGICAL,OPTIONAL :: isPSolve
            END SUBROUTINE pbcg_ilu
         END INTERFACE
         INTERFACE
            SUBROUTINE csr_cg_solver(epsFactor,ncell,neq,ncell_pad,maxmt,maxmt_pad,maxmt_lu0,maxmt_lu1, &
                                     diag,au,ia,ja,ju,                                                  &
                                     diag_lu,alu0,alu1,ia0,ia1,ja0,ja1,                                 &
                                     ap,iap,jap,jaar,iaa,ngroup,nbgroup,                                &
                                     lev_typet,perm_r,                                                  &
                                     source,v1p,izone,isPSolve)
            INTEGER :: ncell,neq,ncell_pad,maxmt,maxmt_pad,maxmt_lu0,maxmt_lu1
            INTEGER :: lev_typet
            INTEGER :: izone
            INTEGER :: ia(ncell+1),ja(maxmt),ju(ncell)
            INTEGER :: ia0(ncell+1),ia1(ncell+1),ja0(maxmt_lu0),ja1(maxmt_lu1)
            INTEGER :: iap(2,ngroup+1),jap(maxmt_pad),jaar(ncell),iaa(2,ngroup+1)
            INTEGER :: ngroup,nbgroup(3,ngroup)
            INTEGER :: perm_r(ncell)
            REAL(8) :: epsFactor
            REAL(8) :: diag(ncell),au(maxmt)
            REAL(8) :: ap(maxmt_pad)
            REAL(8) :: diag_lu(ncell),alu0(maxmt_lu0),alu1(maxmt_lu1)
            REAL(8) :: source(ncell)
            REAL(8) :: v1p(neq)
            LOGICAL, OPTIONAL :: isPSolve
            END SUBROUTINE csr_cg_solver
         END INTERFACE
      END MODULE Zinterface
