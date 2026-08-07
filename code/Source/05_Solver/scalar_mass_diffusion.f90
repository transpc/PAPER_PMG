      SUBROUTINE scalar_mass_diffusion
!     
!     Implement scalar_mass_diffusion 2015.07.29 JHLee (SNU)
!     This routine calculates diffusive fluxes through the cell face
!     between cells i and j.
!     DIFfusive fluxes are discretized using central differences.
!
      USE Zinterface
      USE VOL_DATA     , ONLY: cell
      USE Zvec_param   , ONLY: nf_non,nf_inl
      USE Znum_cell    , ONLY: istart_nf,istart_nbcon_nf, &
                               nf_number_nb,lens,nf_number_id,istart_nfs
      USE Zvec_index   , ONLY: left_nf,right_non,nbcon_nf
      USE Zmass_diff   , ONLY: mdiff_gas,mediff_gas
      USE Zface        , ONLY: laminar
      USE Zb_condition , ONLY: qualab,rhob_gas,alphab_gas
      USE Zvec_geo     , ONLY: f0,f1,sap_nf
      USE Zrv_model    , ONLY: rv_valve      
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i,k
      INTEGER :: ii,kk
      INTEGER :: nv,nf_number,len,istart0,istart,istart2,i0,i1,i2
      REAL(8) :: xn_i,xn_k,hag_i,hag_k
      REAL(8) :: dgi,edgi
      REAL(8) :: xn
!.....Local vector arrays
      REAL(8) :: fluxg_diff_nf(nf_non+nf_inl)
      REAL(8) :: efluxg_diff_nf(nf_non+nf_inl)
!
!.....Build summation info for non,inl
!
      nf_number_nb=1
      nf_number_id(0)=0
      nf_number_id(1)=2
      istart_nfs(0)=0
      istart_nfs(1)=istart_nfs(0)+nf_non
      lens         =istart_nfs(1)+nf_inl
!              
!.....Cells non
!
      nv=0
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      DO i=1,len  
         i0=istart0+i
         i1=istart+i
         ii=left_nf(i1)
         kk=right_non(i)
         xn_i =cell%alphag(ii)*cell%rhog(ii)*cell%mdiff(ii)
         xn_k =cell%alphag(kk)*cell%rhog(kk)*cell%mdiff(kk)
         hag_i=cell%ha(ii)-cell%hg(ii)
         hag_k=cell%ha(kk)-cell%hg(kk)
         dgi =f1(i)*xn_i+f0(i)*xn_k
         edgi=f1(i)*xn_i*hag_i+f0(i)*xn_k*hag_k
         xn =(cell%quala_o(kk)-cell%quala_o(ii))*sap_nf(i1)
         fluxg_diff_nf(i0) =dgi *xn
         efluxg_diff_nf(i0)=edgi*xn
      ENDDO
!
!.....valve model
!       
       IF(rv_valve.eq.1) CALL valve_model_scalar_mass_diffusion(fluxg_diff_nf,efluxg_diff_nf)
!
!.....Cells inl
!
      nv=1
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      istart2=istart_nbcon_nf(nf_number)
      len   =istart_nf(2,nf_number)
      DO i=1,len  
         i0=istart0+i
         i1=istart+i
         i2=istart2+i
         ii=left_nf(i1)
         k=nbcon_nf(i2)
         hag_i=cell%ha(ii)-cell%hg(ii)
         dgi  =alphab_gas(k)*rhob_gas(k)*cell%mdiff(ii)
         edgi =dgi*hag_i
         xn =(qualab(k)-cell%quala_o(ii))*sap_nf(i1)
         fluxg_diff_nf(i0) =dgi *xn
         efluxg_diff_nf(i0)=edgi*xn
      ENDDO
!
      CALL sum_nf(0,-1,                      &
                  fluxg_diff_nf ,mdiff_gas,  &
                  efluxg_diff_nf,mediff_gas)
!
      END SUBROUTINE scalar_mass_diffusion
