!
      SUBROUTINE csr_cg_solver(epsFactor,ncell,neq,ncell_pad,maxmt,maxmt_pad,maxmt_lu0,maxmt_lu1, &
                               diag,au,ia,ja,ju,                                                  &
                               diag_lu,alu0,alu1,ia0,ia1,ja0,ja1,                                 &
                               ap,iap,jap,jaar,iaa,ngroup,nbgroup,                                &
                               lev_typet,perm_r,                                                  &
                               source,v1p,izone,isPSolve)
!
!     This routine select conjugate gradient methods.
!     It is called twice (pressure_solve.f90, calc_solid.f90)
!
!     psolve=1 : conjugate gradient with diagonal precondition
!     psolve=2 : conjugate gradient with ILU precondition
!     psolve=3 : bi-conjugate gradient with diagonal precondition
!     psolve=4 : bi-conjugate gradient with ILU precondition
!
      USE Zinterface
      USE Zbicg , ONLY: psolve,eps_bicg,max_bicg
      USE Zcore , ONLY: myrank,np
!
      IMPLICIT NONE
!
!.....Input
      INTEGER :: ncell,neq,ncell_pad,maxmt,maxmt_pad,maxmt_lu0,maxmt_lu1
      INTEGER :: lev_typet
      INTEGER :: izone
      INTEGER :: ia(ncell+1),ja(maxmt),ju(ncell)
      INTEGER :: ia0(ncell+1),ia1(ncell+1),ja0(maxmt_lu0),ja1(maxmt_lu1)
      INTEGER :: iap(2,ngroup+1),jap(maxmt_pad),iaa(2,ngroup+1)
      INTEGER :: jaar(ncell)
      INTEGER :: ngroup,nbgroup(3,ngroup)
      INTEGER :: perm_r(ncell)
      REAL(8) :: epsFactor
      REAL(8) :: diag(ncell),au(maxmt)
      REAL(8) :: ap(maxmt_pad)
      REAL(8) :: diag_lu(ncell),alu0(maxmt_lu0),alu1(maxmt_lu1)
      REAL(8) :: source(ncell)
      LOGICAL, OPTIONAL :: isPSolve
!.....Output
      REAL(8) :: v1p(neq)
!.....Local variables 
      LOGICAL :: isPSolve_
      REAL(8) :: eps2
!
      if(PRESENT(isPSolve)) then
         isPSolve_ = isPSolve
      else
         isPSolve_ = .false.
      end if
!
      eps2= eps_bicg*epsFactor
!
      SELECT CASE(psolve)
      CASE(1) 
!
!........CG-Diagonal
!      
         CALL pcg_diag(eps2,max_bicg,ncell,neq,maxmt,myrank, &
                       ia,ja,ju,au,source,v1p,np,            &
                       izone)
!         
      CASE(2)
!
!........CG-ILU
!            
         CALL pcg_ilu(eps2,max_bicg,ncell,neq,maxmt,myrank,                        &
                      maxmt_lu0,maxmt_lu1,                                         &
                      ia,ja,ju,au,source,v1p,np,diag_lu,alu0,alu1,ja0,ja1,ia0,ia1, &
                      izone)
!         
      CASE(3)
!
!........BICG-Diagonal
!            
         CALL pbcg_diag(eps2,max_bicg,ncell,neq,maxmt,myrank,  &
                        ia,ja,ju,au,source,v1p,np,             &
                        izone)
!         
      CASE(4)
!
!........BICG-ILU
!            
         CALL pbcg_ilu(eps2,max_bicg,ncell,ncell_pad,maxmt_pad,maxmt_lu0,maxmt_lu1, &
                       diag,                                                        &
                       diag_lu,alu0,alu1,ia0,ia1,ja0,ja1,                           &
                       iap,jap,ap,jaar,iaa,ngroup,nbgroup,                          &
                       lev_typet,perm_r,                                            &
                       neq,source,v1p,izone,isPSolve_)
!         
      END SELECT
!
      END SUBROUTINE csr_cg_solver
