      MODULE Zmars
! 
      IMPLICIT NONE
      SAVE
!
      INTEGER :: n_marsbc,ncupvol,ncell_old(3)
      INTEGER :: n_marsbc_loc                         !mcc-pik-2015-10-21-offdia
      INTEGER,DIMENSION(:),ALLOCATABLE :: cupvol      !mcc-pik-2015-10-21-offdia     
!pik-mcc-mpi-2013-06-26-beg
      INTEGER,DIMENSION(:),ALLOCATABLE :: marsindex
      INTEGER,DIMENSION(:),ALLOCATABLE :: i3invtbl_tmp
      REAL(8) :: time_mars
      REAL(8),DIMENSION(:),ALLOCATABLE :: ppcup,pcup,egcup,elcup,alphagcup, &
                                          rhogcup,rholcup,cboroncup,qualacup
      REAL(8),DIMENSION(:),ALLOCATABLE :: ppcup_tmp,pcup_tmp,egcup_tmp,elcup_tmp,alphagcup_tmp,&
                                          rhogcup_tmp,rholcup_tmp,cboroncup_tmp,qualacup_tmp
      REAL(8),DIMENSION(:),ALLOCATABLE :: tlcup,tgcup,tlcup_tmp,tgcup_tmp
      !The memory for these variables will be allocated in set_cupid_cell_no.f90                             
!pik-mcc-mpi-2013-06-26-end      
      INTEGER,DIMENSION(:),ALLOCATABLE :: i3cupid_loc
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: poiss_mc
!.....Rx trip
      INTEGER :: rx_trip_signal
      REAL(8) :: time_scram
!
!.....CUPID-SG
!.....Heat structure coupling for Steam Generator
      INTEGER,DIMENSION(:),ALLOCATABLE :: nn_marscell
      REAL(8) :: temp_marscell(13),temp_marscell_all(13),temp_marscell_ave(13),tsol_tmp(13)
!   
!.....mass flow rate
      REAL(8),ALLOCATABLE :: mass_nf_mcc(:)
      END MODULE Zmars
