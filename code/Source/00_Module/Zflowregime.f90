      MODULE Zflowregime
! 
      IMPLICIT NONE
      SAVE
!
!     Extremely simplIFied flow regime map
!     Defined by void fraction only.
!     RegimeNo: 1 SPL; 2 BBL; 3 SLC; 4 ANM; 5 MST; 6 SPV
!
      REAL(8),PARAMETER::bbl_slc=0.3D0,slc_anm=0.65D0,anm_mst=0.95D0
      REAL(8) alphag_bc,alphag_cm,gamma_1,gamma_2
      REAL(8) vFgl_1_min,vFgl_2_min,vFgl_3_min
!
      END MODULE Zflowregime
