!
      SUBROUTINE plum_outlet_property_user
!
!     This routine declares flow properties at pressure outlet boundaries
!
      USE Zinterface
      USE VOL_DATA     , ONLY: cell            
      USE Zzone        , ONLY: ncell_fluid
      USE Zparam       , ONLY: ndim
      USE Zvec_param   , ONLY: nf_nonk,nf_non
      USE Znum_cell    , ONLY: istart_nf,                                &
                               right_nb_k,                               &
                               lens,nf_number_nb,nf_number_id,istart_nfs
      USE Zvec_index   , ONLY: left_nf,right_non
      USE Zb_condition , ONLY: alpha_liq_nd,alpha_gas_nd,e_liq_nd,e_gas_nd,quala_nd
      USE Zvector      , ONLY: vl_o,vg_o,vd_o
      USE Zbc_index    , ONLY: npb
      USE Zvec_geo     , ONLY: dji_nf
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i,k
      INTEGER :: ii,kk,ix
      INTEGER :: nv,nf_number,istart0,istart,len,i0,i1
      REAL(8) wsum,wsumr,el,eg,al,ag,quala
      REAL(8) rdji
!.....Local arrays
      REAL(8),DIMENSION(ncell_fluid) :: vel,veg,val,vag,vquala,vwsum
!.....Local vector arrays
      REAL(8),DIMENSION(nf_nonk+nf_non) :: el_nf,eg_nf,al_nf,ag_nf,quala_nf,rdji_nf
!
!.....Build summation info for non
!
      nf_number_nb=0
      nf_number_id(-1)=-1
      nf_number_id(0)=0
      istart_nfs(-1)=0
      istart_nfs(0)=istart_nfs(-1)+nf_nonk
      lens         =istart_nfs(0)+nf_non
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
         IF(npb(ii).gt.0 .and. npb(kk).eq.0) THEN
            rdji=dji_nf(i1)
            el_nf(i0)   =rdji*cell%el(kk)
            eg_nf(i0)   =rdji*cell%eg(kk)
            al_nf(i0)   =rdji*cell%alphal(kk)
            ag_nf(i0)   =rdji*cell%alphag(kk)
            quala_nf(i0)=rdji*cell%quala(kk)
            rdji_nf(i0) =rdji
         ELSE
            el_nf(i0)   =0.d0
            eg_nf(i0)   =0.d0
            al_nf(i0)   =0.d0
            ag_nf(i0)   =0.d0
            quala_nf(i0)=0.d0
            rdji_nf(i0) =0.d0
         ENDIF
      ENDDO
!
      nv=-1
      nf_number=nf_number_id(nv)
      len   =istart_nf(2,nf_number)
      DO i=1,len
         k=right_nb_k(i)
         ii=right_non(k)
         kk=left_nf(k)
         IF(npb(ii).gt.0 .and. npb(kk).eq.0) THEN
            rdji=dji_nf(k)
            el_nf(i)   =rdji*cell%el(kk)
            eg_nf(i)   =rdji*cell%eg(kk)
            al_nf(i)   =rdji*cell%alphal(kk)
            ag_nf(i)   =rdji*cell%alphag(kk)
            quala_nf(i)=rdji*cell%quala(kk)
            rdji_nf(i) =rdji
         ELSE
            el_nf(i)   =0.d0
            eg_nf(i)   =0.d0
            al_nf(i)   =0.d0
            ag_nf(i)   =0.d0
            quala_nf(i)=0.d0
            rdji_nf(i) =0.d0
         ENDIF
      ENDDO
!
      CALL sum_nf(0,0,             &
                  el_nf   ,vel,    &
                  eg_nf   ,veg,    &
                  al_nf   ,val,    &
                  ag_nf   ,vag,    &
                  quala_nf,vquala, &
                  rdji_nf ,vwsum)
!
      DO ix=1,ndim
         DO i=1,ncell_fluid
            IF(npb(i).eq.0) cycle
            vl_o(i,ix)=0.d0
            vg_o(i,ix)=0.d0
            vd_o(i,ix)=0.d0
         ENDDO
      ENDDO
!
      DO i=1,ncell_fluid
         k=npb(i)
         IF(k.eq.0) cycle
         wsum=vwsum(i)
         el=vel(i)
         eg=veg(i)
         al=val(i)
         ag=vag(i)
         quala=vquala(i)
!
!...........Neumann momentum B.C
!...........alpha_gas_nd(1) > 1.0d1 : Press. B.C. + Neumann scalar B.C.
!
         IF(alpha_gas_nd(npb(i)).gt.1.0d1)THEN
            IF(wsum.ne.0.0d0)THEN
               wsumr=1.d0/wsum
               cell%el(i)    =el*wsumr
               cell%eg(i)    =eg*wsumr
               cell%alphal(i)=al*wsumr
               cell%alphag(i)=ag*wsumr
               cell%quala(i) =quala*wsumr
            ELSE
               cell%el(i)    =el
               cell%eg(i)    =eg
               cell%alphal(i)=al
               cell%alphag(i)=ag
               cell%quala(i) =quala
            ENDIF
         ELSE
            cell%el(i)    =e_liq_nd(k)
            cell%eg(i)    =e_gas_nd(k)
            cell%alphal(i)=alpha_liq_nd(k)
            cell%alphag(i)=alpha_gas_nd(k)
            cell%quala(i) =quala_nd(k)
         ENDIF
      ENDDO
!
      END SUBROUTINE plum_outlet_property_user
