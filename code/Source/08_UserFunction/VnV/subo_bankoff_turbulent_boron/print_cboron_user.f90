      SUBROUTINE print_cboron_user(time) 
!
!     Save Output for 3D a graphic processor using Open-GL library
!
      USE Zinterface
      USE VOL_DATA     , ONLY: cell
      USE Zzone        , ONLY: ncell_fluid
      USE Zcore        , ONLY: np,myrank
      USE Zparam       , ONLY: ndim
      USE Zvec_param   , ONLY: nf_inl,nf_out
      USE Znum_cell    , ONLY: istart_nf,istart_nbcon_nf, &
                               nf_number_nb,lens,nf_number_id,istart_nfs
      USE Zvec_index   , ONLY: left_nf,nbcon_nf
      USE Zcoord3      , ONLY: vol   
      USE Zb_condition , ONLY: vb_liq,vb_gas,vb_drp,alphab_gas, &
                               rhob_liq,rhob_gas,eb_liq,eb_gas
      USE Zvector      , ONLY: vl_n,vg_n,vd_n
      USE Zboron       , ONLY: cboronb_liq
      USE Zvec_geo     , ONLY: svp_nf,svp_nf,svp_nf
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i,k
      INTEGER :: ii
      INTEGER :: nf_number,istart,len,istart2,i1,i2
      LOGICAL,SAVE :: initial_p=.true.
      REAL(8) :: time
      REAL(8) :: al,ag,ad
      REAL(8) :: fluxvol_l,fluxvol_d,fluxvol_g
      REAL(8) :: err_boron1,err_boron2,err_boron3
      REAL(8) :: err_mass1,err_mass2,err_mass3
      REAL(8) :: err_erg1,err_erg2,err_erg3
!.....Local arrays
      CHARACTER*25 :: title(25)
      REAL(8) :: tmp(9)
      REAL(8) :: err_boron1v(ncell_fluid),err_mass1v(ncell_fluid),err_erg1v(ncell_fluid)
      REAL(8) :: err_boron2v(ncell_fluid),err_mass2v(ncell_fluid),err_erg2v(ncell_fluid)
!.....Local vector arrays
      REAL(8) :: err_boron_inl(nf_inl),err_mass_inl(nf_inl),err_erg_inl(nf_inl)
      REAL(8) :: err_boron_out(nf_out),err_mass_out(nf_out),err_erg_out(nf_out)
!
      IF(initial_p)THEN
         initial_p = .false.
         IF(myrank.eq.0) OPEN(69,file='VFS12_ref.dat')
         title(1)='A:Time         '  
         title(2)='B:InBoronCon   ' 
         title(3)='C:OutBoronCon  ' 
         title(4)='D:InvBoronCon  ' 
         title(5)='E:InLiquidMass '  
         title(6)='F:OutLiquidMass'  
         title(7)='G:InvLiquidMass' 
         title(8)='H:InLiquidMErg '  
         title(9)='I:OutLiquidErg '  
         title(10)='J:InvLiquidErg ' 
         IF(myrank.eq.0) WRITE(69, 8003)(title(i), i=1,10)
8003  FORMAT(1x,25A20)
      ENDIF
!
      IF(ndim.eq.2) THEN
!
!.....Inlet Flux
!
         nf_number=2
         istart=istart_nf(1,nf_number)
         istart2=istart_nbcon_nf(nf_number)
         len   =istart_nf(2,nf_number)
         DO i=1,len  
            i1=istart+i
            i2=istart2+i
            ii=left_nf(i1)
            k=nbcon_nf(i2)
            fluxvol_l=vb_liq(k,1)*svp_nf(i1,1)+vb_liq(k,2)*svp_nf(i1,2)
            fluxvol_g=vb_gas(k,1)*svp_nf(i1,1)+vb_gas(k,2)*svp_nf(i1,2)
            fluxvol_d=vb_drp(k,1)*svp_nf(i1,1)+vb_drp(k,2)*svp_nf(i1,2)
            err_boron_inl(i)= cboronb_liq(k)*(1.d0-alphab_gas(k))*rhob_liq(k)*fluxvol_l 
            err_mass_inl(i) = (1.d0-alphab_gas(k))               *rhob_liq(k)*fluxvol_l &
                             +      alphab_gas(k)                *rhob_gas(k)*fluxvol_l
            err_erg_inl(i)  = eb_liq(k)*(1.d0-alphab_gas(k))     *rhob_liq(k)*fluxvol_l &
                             +eb_gas(k)*      alphab_gas(k)      *rhob_gas(k)*fluxvol_g
         ENDDO
!
!.....Outlet Flux
!
         nf_number=3
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
         DO i=1,len  
            i1=istart+i
            ii=left_nf(i1)
            fluxvol_l=vl_n(ii,1)*svp_nf(i1,1)+vl_n(ii,2)*svp_nf(i1,2)
            fluxvol_g=vg_n(ii,1)*svp_nf(i1,1)+vg_n(ii,2)*svp_nf(i1,2)
            fluxvol_d=vd_n(ii,1)*svp_nf(i1,1)+vd_n(ii,2)*svp_nf(i1,2)
            err_boron_out(i)= cell%cboron(ii)*cell%alphal(ii)*cell%rhol(ii)*fluxvol_l &
                             +cell%cboron(ii)*cell%alphad(ii)*cell%rhod(ii)*fluxvol_d 
            err_mass_out(i) = cell%alphal(ii)                *cell%rhol(ii)*fluxvol_l &
                             +cell%alphag(ii)                *cell%rhog(ii)*fluxvol_g &
                             +cell%alphad(ii)                *cell%rhod(ii)*fluxvol_d  
            err_erg_out(i)  = cell%el(ii)    *cell%alphal(ii)*cell%rhol(ii)*fluxvol_l &
                             +cell%eg(ii)    *cell%alphag(ii)*cell%rhog(ii)*fluxvol_g &
                             +cell%el(ii)    *cell%alphad(ii)*cell%rhod(ii)*fluxvol_d  
         ENDDO
      ELSE
!
!.....Inlet Flux
!
         nf_number=2
         istart=istart_nf(1,nf_number)
         istart2=istart_nbcon_nf(nf_number)
         len   =istart_nf(2,nf_number)
         DO i=1,len  
            i1=istart+i
            i2=istart2+i
            ii=left_nf(i1)
            k=nbcon_nf(i2)
            fluxvol_l=vb_liq(k,1)*svp_nf(i1,1)+vb_liq(k,2)*svp_nf(i1,2)+vb_liq(k,3)*svp_nf(i1,3)
            fluxvol_g=vb_gas(k,1)*svp_nf(i1,1)+vb_gas(k,2)*svp_nf(i1,2)+vb_gas(k,3)*svp_nf(i1,3)
            fluxvol_d=vb_drp(k,1)*svp_nf(i1,1)+vb_drp(k,2)*svp_nf(i1,2)+vb_drp(k,3)*svp_nf(i1,3)
            err_boron_inl(i)= cboronb_liq(k)*(1.d0-alphab_gas(k))*rhob_liq(k)*fluxvol_l
            err_mass_inl(i) = (1.d0-alphab_gas(k))               *rhob_liq(k)*fluxvol_l &
                             +      alphab_gas(k)                *rhob_gas(k)*fluxvol_l
            err_erg_inl(i)  = eb_liq(k)*(1.d0-alphab_gas(k))     *rhob_liq(k)*fluxvol_l &
                             +eb_gas(k)*      alphab_gas(k)      *rhob_gas(k)*fluxvol_g
         ENDDO
!
!.....Outlet Flux
!
         nf_number=3
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
         DO i=1,len  
            i1=istart+i
            ii=left_nf(i1)
            fluxvol_l=vl_n(ii,1)*svp_nf(i1,1)+vl_n(ii,2)*svp_nf(i1,2)+vl_n(ii,3)*svp_nf(i1,3)
            fluxvol_g=vg_n(ii,1)*svp_nf(i1,1)+vg_n(ii,2)*svp_nf(i1,2)+vg_n(ii,3)*svp_nf(i1,3)
            fluxvol_d=vd_n(ii,1)*svp_nf(i1,1)+vd_n(ii,2)*svp_nf(i1,2)+vd_n(ii,3)*svp_nf(i1,3)
            err_boron_out(i)= cell%cboron(ii)*cell%alphal(ii)*cell%rhol(ii)*fluxvol_l
            err_mass_out(i) = cell%alphal(ii)                *cell%rhol(ii)*fluxvol_l &
                             +cell%alphag(ii)                *cell%rhog(ii)*fluxvol_g
            err_erg_out(i)  = cell%el(ii)    *cell%alphal(ii)*cell%rhol(ii)*fluxvol_l &
                             +cell%eg(ii)    *cell%alphag(ii)*cell%rhog(ii)*fluxvol_g
         ENDDO
      ENDIF
!
!.....Build summation info for inl
!
      nf_number_nb=0
      nf_number_id(0)=2
      istart_nfs(0)=0
      lens         =istart_nfs(0)+nf_inl
!
      CALL sum_nf(0,-1,                      &
                  err_boron_inl,err_boron1v, &
                  err_mass_inl,err_mass1v,   &
                  err_erg_inl,err_erg1v)
!
!.....Build summation info for out
!
      nf_number_nb=0
      nf_number_id(0)=3
      istart_nfs(0)=0
      lens         =istart_nfs(0)+nf_out
      CALL sum_nf(0,-1,                      &
                  err_boron_out,err_boron2v, &
                  err_mass_out,err_mass2v,   &
                  err_erg_out,err_erg2v)
!
      err_boron1=0.0d0
      err_boron2=0.0d0
      err_mass1=0.0d0
      err_mass2=0.0d0
      err_erg1=0.0d0
      err_erg2=0.0d0
      DO i=1,ncell_fluid
         err_boron1=err_boron1+err_boron1v(i)
         err_mass1 =err_mass1 +err_mass1v(i)
         err_erg1  =err_erg1  +err_erg1v(i)
         err_boron2=err_boron2+err_boron2v(i)
         err_mass2 =err_mass2 +err_mass2v(i)
         err_erg2  =err_erg2  +err_erg2v(i)
      ENDDO
!
!.....Inventory
!
      err_boron3=0.0d0
      err_mass3=0.0d0
      err_erg3=0.0d0
      DO i=1,ncell_fluid
         al=cell%alphal(i)*cell%rhol(i)*Vol(i)
         ag=cell%alphag(i)*cell%rhog(i)*Vol(i)
         ad=cell%alphad(i)*cell%rhod(i)*Vol(i) 
         err_boron3=err_boron3 +cell%cboron(i)*al              +cell%cboron(i)*ad
         err_mass3 =err_mass3  +               al+           ag+               ad
         err_erg3  =err_erg3   +cell%el(i)    *al+cell%eg(i)*ag+cell%el(i)    *ad
! bug cell%ed replaced by cell%el
!                            +cell%ed(i)*cell%alphad(i)*cell%rhod(i)*Vol(i) 
      ENDDO
!       
      IF(ABS(err_boron1).le.1.d-99) err_boron1=0.d0      !boron-in
      IF(ABS(err_boron2).le.1.d-99) err_boron2=0.d0      !boron-out
      IF(ABS(err_boron3).le.1.d-99) err_boron3=0.d0      !boron-domain
      IF(ABS(err_mass1) .le.1.d-99) err_mass1 =0.d0 
      IF(ABS(err_mass2) .le.1.d-99) err_mass2 =0.d0 
      IF(ABS(err_mass3) .le.1.d-99) err_mass3 =0.d0 
      IF(ABS(err_erg1)  .le.1.d-99) err_erg1  =0.d0 
      IF(ABS(err_erg2)  .le.1.d-99) err_erg2  =0.d0 
      IF(ABS(err_erg3)  .le.1.d-99) err_erg3  =0.d0
!   
      IF(np.gt.1) THEN
         tmp(1)=err_boron1
         tmp(2)=err_boron2
         tmp(3)=err_boron3
         tmp(4)=err_mass1
         tmp(5)=err_mass2
         tmp(6)=err_mass3
         tmp(7)=err_erg1
         tmp(8)=err_erg2
         tmp(9)=err_erg3
         CALL allreducei_r(tmp,9)
         err_boron1=tmp(1)
         err_boron2=tmp(2)
         err_boron3=tmp(3)
         err_mass1 =tmp(4)
         err_mass2 =tmp(5)
         err_mass3 =tmp(6)
         err_erg1  =tmp(7)
         err_erg2  =tmp(8)
         err_erg3  =tmp(9)
      ENDIF
!                  
      IF(myrank.eq.0) WRITE(69,1000) time,-err_boron1,err_boron2,err_boron3, &
                                          -err_mass1,err_mass2,err_mass3,    &
                                          -err_erg1,err_erg2,err_erg3
1000  FORMAT(1x,10e25.15)
!      
      END SUBROUTINE print_cboron_user    
