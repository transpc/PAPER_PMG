!
      SUBROUTINE IAT_bubsize_change
!            
!     This routine calculates bubble size variation term for IAT.
!     IAT=Interfacial Area Concentration Transport Equation.
!  
      USE Zinterface
      USE VOL_DATA              
      USE Zmpi         , ONLY: ncell_fp
      USE Zzone        , ONLY: ncell_fluid
      USE Zvec_param   , ONLY: nf_non,nf_mcc,nf_inl,nf_out,nf_flux
      USE Znum_cell    , ONLY: istart_nf,istart_nbcon_nf, &
                               nf_number_nb,lens,nf_number_id,istart_nfs
      USE Zvec_index   , ONLY: left_nf,right_non,nbcon_nf
      USE Zconst2      , ONLY: dt
      USE Ziat         , ONLY: ia_old,iat_size
      USE Zqvol        , ONLY: gamma
      USE Zb_condition , ONLY: rhob_gas
      USE Zvec_major   , ONLY: flux_g_nf
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i,k
      INTEGER :: ii,kk
      INTEGER :: nv,nf_number,istart,isize,istart2,i1,i2
      REAL(8) :: rhog1,rhog2
      REAL(8) :: dsmi,coeff,ia_condens,ia_pdrop 
      REAL(8) :: ia_mass_fluxs,ia_fluxs
!.....Local arrays
      REAL(8) :: ia_cond_pdrop(ncell_fp)
      REAL(8) :: ia_mass_flux(ncell_fluid),ia_flux(ncell_fluid)
!.....Local vector arrays
      REAL(8) :: flux_conv_nf(nf_flux)
!
!.....Build summation info for non,mcc,inl,out
!
      nf_number_nb=3
      nf_number_id(0)=0
      nf_number_id(1)=1
      nf_number_id(2)=2
      nf_number_id(3)=3
      istart_nfs(0)=0
      istart_nfs(1)=istart_nfs(0)+nf_non
      istart_nfs(2)=istart_nfs(1)+nf_mcc
      istart_nfs(3)=istart_nfs(2)+nf_inl
      lens         =istart_nfs(3)+nf_out
!
!.....1st and 2nd terms of size variation term
!
      DO i=1,ncell_fluid
         dsmi=cell%d1(i)
         IF(cell%alphag(i).gt.1.d-3)THEN
            coeff=2.0d0/3.d0*ia_old(i)/(cell%alphag(i)*cell%rhog(i))
         ELSE
            coeff=0.0d0
         ENDIF
         IF(gamma(i).lt.0.d0)ia_condens=coeff*gamma(i) 
         ia_condens=coeff*gamma(i)
         ia_pdrop=coeff*cell%alphag(i)*(cell%rhog(i)-cell%rhog_o(i))/dt
         ia_cond_pdrop(i)=ia_condens-ia_pdrop
      ENDDO
!
!.....3rd terms of size variation term
!
!.....Cells non
!
      nv=0
      nf_number=nf_number_id(nv)
      istart=istart_nf(1,nf_number)
      isize =istart_nf(2,nf_number)
      DO i=1,isize
         i1=istart+i
         ii=left_nf(i1)
         kk=right_non(i)
         rhog1=cell%rhog(ii)
         rhog2=cell%rhog(kk)
         flux_conv_nf(i1)=MIN(flux_g_nf(i1),0.d0)*rhog2+MAX(flux_g_nf(i1),0.d0)*rhog1
      ENDDO
!
!.....Celles mcc
!
      nv=1
      nf_number=nf_number_id(nv)
      istart=istart_nf(1,nf_number)
      isize =istart_nf(2,nf_number)
      DO i=1,isize
         i1=istart+i
         ii=left_nf(i1)
         rhog1=cell%rhog(ii)
         flux_conv_nf(i1)=flux_g_nf(i1)*rhog1
      ENDDO
!
!.....Cells inl
!
      nv=2
      nf_number=nf_number_id(nv)
      istart=istart_nf(1,nf_number)
      isize =istart_nf(2,nf_number)
      istart2=istart_nbcon_nf(nf_number)
      DO i=1,isize
         i1=istart+i
         i2=istart2+i
         ii=left_nf(i1)
         k=nbcon_nf(i2)
         rhog1=cell%rhog(ii)
         rhog2=rhob_gas(k)
         flux_conv_nf(i1)=MIN(flux_g_nf(i1),0.d0)*rhog2+MAX(flux_g_nf(i1),0.d0)*rhog1
      ENDDO
!
!.....Cells out
!
      nv=3
      nf_number=nf_number_id(nv)
      istart=istart_nf(1,nf_number)
      isize =istart_nf(2,nf_number)
      DO i=1,isize
         i1=istart+i
         ii=left_nf(i1)
         rhog1=cell%rhog(ii)
         flux_conv_nf(i1)=flux_g_nf(i1)*rhog1
      ENDDO
!
      CALL sum_nf(0,-1,                      &
                  flux_conv_nf,ia_mass_flux, &
                  flux_g_nf   ,ia_flux)
!      
!.....Sum of 1st,2nd,3rd,and 4th terms of size variation term
!
      DO i=1,ncell_fluid
         ia_mass_fluxs=ia_mass_flux(i)
         ia_fluxs=ia_flux(i)
         iat_size(i)=ia_cond_pdrop(i)-2.0d0/3.0d0*ia_old(i)*(1.d0/cell%rhog(i)*ia_mass_fluxs-1.d0*ia_fluxs)
      ENDDO
!
      END SUBROUTINE IAT_bubsize_change
