      MODULE Zparam
!
      IMPLICIT NONE
      SAVE
!
      INTEGER ndim,nn,ns,mesh_openfoam,mesh_binary
!
!.....Parameter for general usage
!
      INTEGER ::nb_max
      INTEGER ::nin_max
      INTEGER ::nb_sym
      INTEGER ::nb_mars
!
!     REAL(8),PARAMETER::pi=3.14159265358979323846d0
      REAL(8),PARAMETER::pi     =z'400921FB54442D18'
      REAL(8),PARAMETER::pio2   =z'3FF921FB54442D18'
      REAL(8),PARAMETER::pi2    =z'401921FB54442D18'
      REAL(8),PARAMETER::pi4    =z'402921FB54442D18'
      REAL(8),PARAMETER::sqrt_pi=z'3FFC5BF891B4EF6A'
      CHARACTER(50) outfilename
!
! Parameter for Turbulence
      REAL(8) cmu,cappa,coeff_B,prt,ced1,ced2
      REAL(8) ke_cff,dp_cff,clog,ke_small,c_td
      REAL(8) ced1_RNG,ced2_RNG,RNG_cff                ! for RNG k-e model
      REAL(8) ced2_Real,ke_cff_Real,dp_cff_Real        ! for Realizable k-e model      
      REAL(8) RNG_cffr,ke_cff_Realr,ke_cffr,dp_cff_Realr,dp_cffr
      REAL(8) prtr
      PARAMETER (cmu=0.09d0,cappa=0.41d0,coeff_B=5.5d0)
      PARAMETER (prt=0.9d0,ced1=1.44d0,ced2=1.92d0)    !water, standard k-e
      PARAMETER (ke_cff=1.0d0,dp_cff=1.3d0)            !standard k-e
      PARAMETER (ced1_RNG=1.42d0,ced2_RNG=1.68d0,RNG_cff=0.7194d0)    !RNG k-e
      PARAMETER (ced2_Real=1.9d0,ke_cff_Real=1.0d0,dp_cff_Real=1.2d0) !Realizable k-e      
      PARAMETER (clog=9.8d0,ke_small=1.e-15)
      PARAMETER (c_td=0.1d0)
      PARAMETER (RNG_cffr=1.d0/RNG_cff)
      PARAMETER (ke_cff_Realr=1.d0/ke_cff_Real,ke_cffr=1.d0/ke_cff)
      PARAMETER (dp_cff_Realr=1.d0/dp_cff_Real,dp_cffr=1.d0/dp_cff)
      PARAMETER (prtr=1.d0/prt)
!
!.....Parameter for IAT
!
      REAL(8) kc1,kc2,kc3,kb1,kb2
      REAL(8) weber_cr,alpha_max
      REAL(8) Kfactor_iat
      PARAMETER (kc1=2.86d0,kc2=1.922d0,kc3=1.017d0)
      PARAMETER (kb1=1.6d0,kb2=0.42d0)
      PARAMETER (weber_cr=1.24d0,alpha_max=0.52d0)
      PARAMETER (Kfactor_iat=1.d0)
!
      END MODULE Zparam
