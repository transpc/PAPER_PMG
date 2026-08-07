!
      SUBROUTINE check_mass
!
!     This routine calculate mass error and total mass of the system
!
      USE Zinterface
      USE VOL_DATA     , ONLY: cell
      USE Zzone        , ONLY: ncell_fluid
      USE Zvec_param   , ONLY: nf_inl,nf_out,nf_mcc
      USE Zvec_index   , ONLY: left_nf,nbcon_nf
      USE Znum_cell    , ONLY: istart_nf,istart_nbcon_nf,                &
                               nf_number_nb,lens,nf_number_id,istart_nfs
      USE Zb_condition , ONLY: alphab_liq,alphab_drp,alphab_gas,rhob_liq,rhob_drp,rhob_gas
      USE Zcoord3      , ONLY: volp
      USE Zmass_conv   , ONLY: tot_mass
      USE Ztimecon     , ONLY: error_mass
      USE c3com_cupid  , ONLY: i3invtbl,mcdirect,c3dpv,mcgdirect
      USE Zvec_major   , ONLY: flux_l_nf,flux_g_nf,flux_d_nf
      USE Zmars        , ONLY: mass_nf_mcc
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i,idx,k
      INTEGER :: ii
      INTEGER :: nv,nf_number,len,istart0,istart,istart2,i0,i1,i2
      REAL(8) :: ardi,argi,arli
!.....Local arrays
      REAL(8) :: error_massv(ncell_fluid)
!.....Local vectorr array
      REAL(8) :: error_mass_nf(nf_mcc+nf_inl+nf_out)
!
!.....Build summation info for mcc,inl,out
!
      nf_number_nb=2
      nf_number_id(0)=1
      nf_number_id(1)=2
      nf_number_id(2)=3
      istart_nfs(0)=0
      istart_nfs(1)=istart_nfs(0)+nf_mcc
      istart_nfs(2)=istart_nfs(1)+nf_inl
      lens         =istart_nfs(2)+nf_out
!
!.....Cells mcc
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
         idx=i3invtbl(i)
!
!........Set interface flow property
!         mcdirect(idx)<0 : MARS  --> CUPID
!         mcdirect(idx)>0 : CUPID --> MARS
!
          IF(mcdirect(idx).lt.0)THEN
             arli=c3dpv(idx,5)
             ardi=c3dpv(idx,4)
          ELSE
             arli=cell%alphal(i)*cell%rhol(i)
             ardi=cell%alphad(i)*cell%rhol(i)
          ENDIF
             IF(mcgdirect(idx).lt.0)THEN
                argi=c3dpv(idx,3)
             ELSE
                argi=cell%alphag(ii)*cell%rhog(ii)
             ENDIF
             error_mass_nf(i0)=arli*flux_l_nf(i1)+argi*flux_g_nf(i1)+ardi*flux_d_nf(i1)
             mass_nf_mcc(i)=error_mass_nf(i0)
      ENDDO
!
!.....Cells inl
!
      nv=1
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      istart2=istart_nbcon_nf(nf_number)
      DO i=1,len  
         i0=istart0+i
         i1=istart+i
         i2=istart2+i
         ii=left_nf(i1)
         k=nbcon_nf(i2)
         error_mass_nf(i0)= alphab_liq(k)*rhob_liq(k)*flux_l_nf(i1)  &
                           +alphab_gas(k)*rhob_gas(k)*flux_g_nf(i1)  &
                           +alphab_drp(k)*rhob_drp(k)*flux_d_nf(i1)
      ENDDO
!
!.....Cells out
!
      nv=2
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      DO i=1,len  
         i0=istart0+i
         i1=istart+i
         ii=left_nf(i1)
         error_mass_nf(i0)= cell%alphal(ii)*cell%rhol(ii)*flux_l_nf(i1) &
                           +cell%alphag(ii)*cell%rhog(ii)*flux_g_nf(i1) &
                           +cell%alphad(ii)*cell%rhod(ii)*flux_d_nf(i1)
      ENDDO
!
      CALL sum_nf(0,-1,                      &
                  error_mass_nf,error_massv)
!
      error_mass=0.0d0
      DO i=1,ncell_fluid
         error_mass=error_mass+error_massv(i)
      ENDDO
!
!........Calculate total mass
!
      tot_mass=0.d0
      DO i=1,ncell_fluid
         tot_mass=tot_mass+volp(i)*( cell%alphal(i)*cell%rhol(i)   &
                                    +cell%alphag(i)*cell%rhog(i)   &
                                    +cell%alphad(i)*cell%rhod(i) )
      ENDDO
!
!.....Convert the error mass in %
!
      IF(tot_mass.gt.0.0d0)THEN
         error_mass=error_mass/tot_mass*100.0d0
      ELSE
         error_mass=0.0d0
      ENDIF
!
      END SUBROUTINE check_mass
