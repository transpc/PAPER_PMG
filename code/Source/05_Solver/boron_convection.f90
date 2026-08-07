!
      SUBROUTINE boron_convection
!
!     This routine calculates boron convective fluxes through the cell face
!
      USE Zinterface
      USE VOL_DATA     , ONLY: cell                
      USE Zmpi         , ONLY: ncell_fp
      USE Zzone        , ONLY: ncell_fluid
      USE Zcore        , ONLY: np
      USE Zb_condition , ONLY: alphab_liq,rhob_liq
      USE Zare         , ONLY: ar_liq
      USE Zboron       , ONLY: cboronb_liq
      USE Zvec_param   , ONLY: nf_flux
      USE Znum_cell    , ONLY: istart_nf,istart_nbcon_nf
      USE Zvec_index   , ONLY: left_nf,right_non,nbcon_nf
      USE Zvec_major   , ONLY: flux_l_nf,lbor_conv_nf
      USE Z2nd_order   , ONLY: boron_conv_2nd
      USE c3com_cupid  , ONLY: i3invtbl,mcdirect,c3dpv
                               
!
      IMPLICIT NONE
!.....Local variables   
      LOGICAL,SAVE :: initial=.true.
      INTEGER :: i,k
      INTEGER :: ii,kk,mm
      INTEGER :: nf_number,istart,isize,istart2,i1,i2,len,idx
      REAL(8) :: arc1_liq,arc2_liq
      REAL(8) :: lbor_conv_up
!.....Local arrays
      REAL(8) :: arc_liq(ncell_fp)
!
!
!.....Save convective variables
!
      DO i=1,ncell_fluid
         arc_liq(i)=ar_liq(i)*cell%cboron(i)
      ENDDO
!
!.....Communicate convective variables
!
      IF(np.gt.1) CALL communicate_1d(arc_liq)
!
      IF(initial)THEN
         ALLOCATE(lbor_conv_nf(nf_flux))
         initial=.FALSE.
      ENDIF
!
!.....Calculate convective flux through cell face
!
!
!.....Computing cells
!
      nf_number=0
      istart=istart_nf(1,nf_number)
      isize =istart_nf(2,nf_number)
      DO i=1,isize
         i1=istart+i
         ii=left_nf(i1)
         kk=right_non(i)
!         
!...........Apply first order upwind convection
!
         IF    (flux_l_nf(i1).gt.0.d0)THEN
            mm=ii
         ELSE
            mm=kk
         ENDIF
         lbor_conv_nf(i1)=arc_liq(mm)
!      
      ENDDO
!
!.....Inlet
!
      nf_number=2
      istart=istart_nf(1,nf_number)
      isize =istart_nf(2,nf_number)
      istart2=istart_nbcon_nf(nf_number)
      DO i=1,isize
         i1=istart+i
         i2=istart2+i
         ii=left_nf(i1)
         k=nbcon_nf(i2)
!         
         arc1_liq=arc_liq(ii)
!
         arc2_liq=alphab_liq(k)*rhob_liq(k)*cboronb_liq(k)
!
         IF    (flux_l_nf(i1).gt.0.d0)THEN
            lbor_conv_up=arc1_liq
         ELSEIF(flux_l_nf(i1).lt.0.d0)THEN
            lbor_conv_up=arc2_liq
         ELSE
            lbor_conv_up=(arc1_liq+arc2_liq)*0.5d0
         ENDIF
!
         lbor_conv_nf(i1)=lbor_conv_up
!         
      ENDDO
!
!.....Outlet
!
      nf_number=3
      istart=istart_nf(1,nf_number)
      isize =istart_nf(2,nf_number)
      DO i=1,isize
         i1=istart+i
         ii=left_nf(i1)
!         
         mm=ii
!
         lbor_conv_nf(i1)=arc_liq(mm)
!
      ENDDO
!
!.....MARS interface
!
      nf_number=1
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      DO i=1,len  
         i1=istart+i
         ii=left_nf(i1)
         arc1_liq=arc_liq(ii)
!         
         idx=i3invtbl(i)
         IF(mcdirect(idx).lt.0)THEN !from mars to cupid
            arc2_liq=c3dpv(idx,8)
         ELSE
            arc2_liq=arc_liq(ii)
         ENDIF
!
         IF    (flux_l_nf(i1).gt.0.d0)THEN
            lbor_conv_up=arc1_liq
         ELSEIF(flux_l_nf(i1).lt.0.d0)THEN
            lbor_conv_up=arc2_liq
         ELSE
            lbor_conv_up=(arc1_liq+arc2_liq)*0.5d0
         ENDIF
         lbor_conv_nf(i1)=lbor_conv_up
      ENDDO         
!
!.....symmetric, wall condition : conv=0               
!
!
      IF(boron_conv_2nd.gt.0) CALL boron_2nd_conv
!		
      END SUBROUTINE boron_convection
