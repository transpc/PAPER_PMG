! ====================================================================!
!  module for MPI
!---------------------------------------------------------------------+
      module MD_MPI
      implicit none
! ---
      integer nnbd,nintf,nintr,nprcs,myrank,myrankt
      integer,dimension(:),allocatable :: spt,rpt,sintf,rintf,nbdom
      integer nnsend,nnrecv
      
!---
      save
    end module
	
! ====================================================================!
!  module for MPI-MG
!---------------------------------------------------------------------+
      module MD_MPI_MG
      implicit none
! ---
      integer nnbd1,nintf1,nnbd2,nintf2
      integer nmaxgl,nnzr1
      integer,dimension(:),allocatable :: spt1,rpt1,sintf1,rintf1,nbdom1
      integer,dimension(:),allocatable :: siaf,riaf
      
! NEW
      INTEGER(4),DIMENSION(:),ALLOCATABLE:: iintf,inodegl
      INTEGER(4) nintfs,nnbds
      INTEGER(4),DIMENSION(:),ALLOCATABLE:: spts,rpts,sintfs,rintfs,nbdoms
      INTEGER(4),DIMENSION(:),ALLOCATABLE:: inbdc
      INTEGER(4),DIMENSION(:,:),ALLOCATABLE:: ibdomc,isptc,irptc,isintfc,irintfc,isiac,iriac
      INTEGER(4) nnsend_m,nnrecv_m
      INTEGER (4) icommu,iGS,nGS,iallocate_c
 
!---
      save
    end module
    
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
!---------------------------------------------------------------------+
      module MD_MPI_ARP
      implicit none
! ---
      integer(4) nnbdA,nnbdR,nnbdP
      integer(4),dimension(:),allocatable :: sptA,rptA,sintfA,rintfA,nbdomA
      integer(4),dimension(:),allocatable :: sptR,rptR,sintfR,rintfR,nbdomR 
      integer(4),dimension(:),allocatable :: sptP,rptP,sintfP,rintfP,nbdomP
      integer(4) nnsendA,nnrecvA
      integer(4) nnsendR,nnrecvR
      integer(4) nnsendP,nnrecvP
!       
      INTEGER(4) nnbdsA,nnbdsR,nnbdsP
      INTEGER(4),DIMENSION(:),ALLOCATABLE:: sptsA,rptsA,sintfsA,rintfsA,nbdomsA
      INTEGER(4),DIMENSION(:),ALLOCATABLE:: sptsR,rptsR,sintfsR,rintfsR,nbdomsR
      INTEGER(4),DIMENSION(:),ALLOCATABLE:: sptsP,rptsP,sintfsP,rintfsP,nbdomsP
      INTEGER(4),DIMENSION(:),ALLOCATABLE:: inbdcA,inbdcR,inbdcP
      INTEGER(4),DIMENSION(:,:),ALLOCATABLE:: ibdomcA,isptcA,irptcA,isintfcA,irintfcA
      INTEGER(4),DIMENSION(:,:),ALLOCATABLE:: ibdomcR,isptcR,irptcR,isintfcR,irintfcR
      INTEGER(4),DIMENSION(:,:),ALLOCATABLE:: ibdomcP,isptcP,irptcP,isintfcP,irintfcP
      INTEGER(4) nnsend_mA,nnrecv_mA
      INTEGER(4) nnsend_mR,nnrecv_mR
      INTEGER(4) nnsend_mP,nnrecv_mP
      
! temporary for serial
      
      INTEGER(4),DIMENSION(:),ALLOCATABLE:: nnbdomA,nnbdomR,nnbdomP,cext_tmp
      INTEGER(4),DIMENSION(:,:),ALLOCATABLE:: inbdomA,riA,siA,rintA,sintA
      INTEGER(4),DIMENSION(:,:),ALLOCATABLE:: inbdomR,riR,siR,rintR,sintR
      INTEGER(4),DIMENSION(:,:),ALLOCATABLE:: inbdomP,riP,siP,rintP,sintP
 
!---
      save
    end module
	
    ! = = = = = = = = = = = = = = = = = = = = =     
      module MD_OpenMP
      implicit none
! ---
      integer :: nthre
      
!---
      save
    end module
	
     
