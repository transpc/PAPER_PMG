!
      SUBROUTINE scalar_matrix
!
!     This routine determines the elements of 6x6 scalar matrix
!
      USE VOL_DATA      , ONLY: cell
      USE Zmpi          , ONLY: ia_nrhs
      USE Zzone         , ONLY: ncell_fluid
      USE Zparam        , ONLY: ns
      USE Zvec_param    , ONLY: nf_nonk,nf_flux1
      USE Znum_cell     , ONLY: istart_nf, &
                                iptr_nb_k,ia_nb,right_nb_k, &
                                ia_nb,istart_nb1,icell_nb
      USE Zvec_index    , ONLY: left_nf,right_non
      USE Zbc_index     , ONLY: npb
      USE Zporous       , ONLY: tm_mas_l,tm_mas_g,tm_eng_l,tm_eng_g,vd_mas_l,vd_mas_g,vd_eng_l,vd_eng_g 
      USE Zconst2       , ONLY: dt
      USE Zcoord3       , ONLY: volpr
      USE Zdhda         , ONLY: dHldtl,dHgdtg
      USE Zenergy_diff  , ONLY: ediff_gas,ediff_liq
      USE Zmass_diff    , ONLY: mdiff_gas,mediff_gas      
      USE Zpress        , ONLY: p
      USE Zqvol         , ONLY: gamma,gamma_wall,h_ig,h_il,h_gf,qvol_gas,qvol_liq, &
                                qporous_gas,qporous_liq
      USE Zscalar_coeff , ONLY: sb,l_th_equil,alphag_min,alphal_min
      USE Zmodel        , ONLY: rad_source,rad_model
      USE Zscalar_coeff , ONLY: sfg_nf,sfl_nf,sfd_nf,           &
                                sfg_non_k,sfl_non_k,sfd_non_k,  &
                                sfg6_nf,sfl6_nf,sfd6_nf,        &
                                sfg6_non_k,sfl6_non_k,sfd6_non_k
      USE Zvec_major    , ONLY: liq_conv_nf,vap_conv_nf,drp_conv_nf,  &
                                ecnvc_l_nf,ecnvc_g_nf,ecnvc_d_nf,     &
                                void_conv_nf,quala_conv_nf
      !OPR1000 rod-scale (mixing vane)
      USE Zporous       , ONLY: mixing_vane_l
      USE Ztimecon      , ONLY: time
      USE Zporous  , ONLY: l_subchannel,l_mixing_vane
!
      IMPLICIT NONE
!
      INTEGER, PARAMETER :: nblk=64
!.....Local variables
      INTEGER :: i,j,k,j1
      INTEGER :: ii,kk,is,ie,nb,i3,ip
      INTEGER :: n_rhs
      INTEGER :: nf_number,istart,len,isize,istart1,i0,i1
      REAL(8) :: hi_gas,hi_liq,del_hi,hi_gas_w,hi_liq_w
      REAL(8) :: ag_eg,dtdhfg
      REAL(8) :: c_ig,c_if,c_if_hg,c_ig_hf
      REAL(8) :: c_ig_dtvde,c_if_dtlde,c_ig_dtvde_hf,c_if_dtlde_hg
      REAL(8) :: del_dtvdp,del_dtldp,del_tsg,del_tsl
      REAL(8) :: al,ad
      REAL(8) :: c_ht,c_pdt,dtvol
      REAL(8) :: dsesd,yeta1
      REAL(8) :: PsP,c_gf,c_ig_dtsde,c_if_dtsde,c_gf_dtvde,c_gf_dtlde
      REAL(8) :: c_ig_dtsde_hf,c_if_dtsde_hg
      REAL(8) :: del_tgl
      REAL(8) :: cig,cif,cgf,cig_hf,cif_hg
      REAL(8) :: a,b,c
!.....Local arrays
      INTEGER :: ip1(0:2),ip2(0:2)
      REAL(8) :: sm(nblk,6,6)
      REAL(8) :: bm(1+3*ns,6,nblk)
      REAL(8) :: bm1(nblk,6)
!.....Local vector arrays
      REAL(8) :: fg_non_k(nf_nonk)
      REAL(8) :: fg_nf(nf_flux1),fl_nf(nf_flux1),fd_nf(nf_flux1),fq_nf(nf_flux1)
!
!#include '../00_Module/c_Solver/avx.h'
!
!      IF(udfl_gamma_linear)THEN
!         alphadiffr=1.d0/(alphag_cm-alphag_bc)
!      ENDIF
!
!.....Scalar update for Mixing Vane
!
      IF(l_mixing_vane .and. time.ge.0.1d0) then
         CALL udfn_subchannel_scalar_conv_mixingvane
      ENDIF      
!
      nf_number=-1
      istart=istart_nf(1,nf_number)
      isize =istart_nf(2,nf_number)
      DO i=1,isize
         k=right_nb_k(i)
         kk=right_non(k)
         fg_non_k(i)=-(ecnvc_g_nf(k)+void_conv_nf(k)*p(kk))
      ENDDO
!
      DO nf_number=0,2
         istart=istart_nf(1,nf_number)
         isize =istart_nf(2,nf_number)
         DO i=1,isize
            i1=istart+i
            ii=left_nf(i1)
            fg_nf(i1)=-(ecnvc_g_nf(i1)+void_conv_nf(i1)*p(ii))
            fl_nf(i1)=-ecnvc_l_nf(i1)
            fd_nf(i1)=-ecnvc_d_nf(i1)
            fq_nf(i1)=-quala_conv_nf(i1)
         ENDDO
      ENDDO
!     
      ip1(0)=1
      ip1(1)=1
      ip1(2)=1
      ip2(0)=1
      ip2(1)=1
      ip2(2)=1
      DO is=1,ncell_fluid,nblk
         ie=min(is+nblk-1,ncell_fluid)
         nb=ie-is+1
         DO i=is,ie
            IF(npb(i).ne.0) CYCLE
            i3=i-is+1
            n_rhs=ia_nrhs(i)
            dtvol=dt*volpr(i)
            DO j=1,3*n_rhs
               bm(1+j,1,i3)=0.d0
               bm(1+j,2,i3)=0.d0
               bm(1+j,3,i3)=0.d0
               bm(1+j,4,i3)=0.d0
               bm(1+j,5,i3)=0.d0
               bm(1+j,6,i3)=0.d0
            ENDDO
            ii=iptr_nb_k(i)
            j=0
            IF(ii.gt.0) THEN
!DIR$ IVDEP
               DO ip=ia_nb(ii),ia_nb(ii+1)-1
                  k=right_nb_k(ip)
                  j=j+1
                  bm(1+        j,1,i3)=fg_non_k(ip)*dtvol
                  bm(1+  n_rhs+j,2,i3)=fl_nf(k)*dtvol
                  bm(1+2*n_rhs+j,2,i3)=fd_nf(k)*dtvol
                  bm(1+        j,3,i3)=-vap_conv_nf(k)*dtvol
                  bm(1+2*n_rhs+j,4,i3)=-drp_conv_nf(k)*dtvol
                  bm(1+  n_rhs+j,5,i3)=-liq_conv_nf(k)*dtvol
                  bm(1+        j,6,i3)=fq_nf(k)*dtvol
               ENDDO
            ENDIF
            DO nf_number=0,2
               ii=ip1(nf_number)
               istart1=istart_nb1(1,nf_number)
               len    =istart_nb1(2,nf_number)
               IF(ii.gt.len) cycle
               i1=istart1+ii
!
!..............Find pointer ii to access _nb that correspond to i cell
!..............and save for next i search
!
300            CONTINUE
               IF(icell_nb(i1).lt.i) THEN
                  ii=ii+1
                  IF(ii.gt.len) goto 310
                  i1=istart1+ii
                  ip1(nf_number)=ii
                  GOTO 300
               ENDIF
               IF(icell_nb(i1).eq.i) THEN
                  istart =istart_nf(1,nf_number)
                  istart1=istart_nb1(1,nf_number)
                  i0=istart1+ii
!DIR$ IVDEP
                  DO ip=ia_nb(i0),ia_nb(i0+1)-1
                     i1=istart+ip
                     j=j+1
                     bm(1+        j,1,i3)=fg_nf(i1)*dtvol
                     bm(1+  n_rhs+j,2,i3)=fl_nf(i1)*dtvol
                     bm(1+2*n_rhs+j,2,i3)=fd_nf(i1)*dtvol
                     bm(1+        j,3,i3)=-vap_conv_nf(i1)*dtvol
                     bm(1+2*n_rhs+j,4,i3)=-drp_conv_nf(i1)*dtvol
                     bm(1+  n_rhs+j,5,i3)=-liq_conv_nf(i1)*dtvol
                     bm(1+        j,6,i3)=fq_nf(i1)*dtvol
                  ENDDO
                  ip1(nf_number)=ip1(nf_number)+1
               ENDIF
310            CONTINUE
            ENDDO
         ENDDO
!      
!DIR$ VECTOR ALIGNED
         DO i=is,ie
            IF(npb(i).ne.0) CYCLE
            i3=i-is+1
            IF(gamma(i).ge.0.d0)THEN
               hi_gas=cell%hgsat(i)
               hi_liq=cell%hl(i)
            ELSE
               hi_gas=cell%hg(i)
               hi_liq=cell%hlsat(i)
            ENDIF
            IF(gamma_wall(i).ge.0.d0)THEN
               hi_gas_w=cell%hgsat(i)
               hi_liq_w=cell%hl(i)
            ELSE
               hi_gas_w=cell%hg(i)
               hi_liq_w=cell%hlsat(i)
            ENDIF            
!
            PsP=cell%pps(i)/p(i)
            ag_eg=cell%alphag(i)*cell%eg(i)
            del_hi=hi_gas-hi_liq
            dtdhfg=dt/del_hi
!
            cig=dtdhfg*PsP
            cif=dtdhfg
            cgf=dt*(1.d0-PsP)
            c_ig=cig*H_ig(i)
            c_if=cif*H_il(i)
            c_gf=cgf*H_gf(i)
!
            c_ig_dtvde=c_ig*cell%dtgde(i)
            c_ig_dtsde=c_ig*cell%dtsde(i)
            c_if_dtlde=c_if*cell%dtlde(i)
            c_if_dtsde=c_if*cell%dtsde(i)
            c_gf_dtvde=c_gf*cell%dtgde(i)
            c_gf_dtlde=c_gf*cell%dtlde(i)
!
            c_ig_dtvde_hf=c_ig_dtvde*hi_liq
            c_ig_dtsde_hf=c_ig_dtsde*hi_liq
            c_if_dtlde_hg=c_if_dtlde*hi_gas
            c_if_dtsde_hg=c_if_dtsde*hi_gas
!
            del_dtvdp=cell%dtsdp(i)-cell%dtgdp(i)
            del_dtldp=cell%dtsdp(i)-cell%dtldp(i)
!
            c_ig_hf=c_ig*hi_liq
            c_if_hg=c_if*hi_gas
            c_pdt=c_ig_hf*del_dtvdp+c_if_hg*del_dtldp
!
!...........Vapor energy equation
!
            a=c_ig_dtsde_hf-c_ig_dtvde_hf+c_if_dtsde_hg+c_gf_dtvde
            b=c_ig_hf*(cell%dtsdx(i)-cell%dtgdx(i))+c_if_hg*cell%dtsdx(i)+c_gf*cell%dtgdx(i)
            c=c_pdt-c_gf*(cell%dtgdp(i)+cell%dtldp(i))
            sm(i3,1,1)=ag_eg*cell%drhogde(i)+cell%alphag(i)*cell%rhog(i)+a
            sm(i3,1,2)=-(c_if_dtlde_hg+c_gf_dtlde)
            sm(i3,1,3)=cell%rhog(i)*cell%eg(i)+p(i)
            sm(i3,1,4)=0.d0
            sm(i3,1,5)= ag_eg*cell%drhogdx(i)+b
            sm(i3,1,6)= ag_eg*cell%drhogdp(i)+c
!
!...........Liquid energy equation
!
!
            al=cell%alphal(i)
            ad=cell%alphad(i)
            sm(i3,2,1)=-a
            sm(i3,2,2)=(al+ad)*cell%el(i)*cell%drholde(i)+(al+ad)*cell%rhol(i)+c_if_dtlde_hg+c_gf_dtlde
            sm(i3,2,3)=-cell%rhol(i)*cell%el(i)-p(i)
            sm(i3,2,4)=0.d0
            sm(i3,2,5)=-b
            sm(i3,2,6)=(al+ad)*cell%el(i)*cell%drholdp(i)-c
!...........Vapor continuity equation
!
            sm(i3,3,1)=cell%alphag(i)*cell%drhogde(i)+c_ig_dtsde-c_ig_dtvde+c_if_dtsde
            sm(i3,3,2)=-c_if_dtlde
            sm(i3,3,3)=cell%rhog(i)
            sm(i3,3,4)=0.d0
            sm(i3,3,5)=cell%alphag(i)*cell%drhogdx(i)+c_ig*(cell%dtsdx(i)-cell%dtgdx(i))+c_if*cell%dtsdx(i)
            sm(i3,3,6)=cell%alphag(i)*cell%drhogdp(i)+c_ig*del_dtvdp+c_if*del_dtldp
!
            dtvol=dt*volpr(i)
            del_tsg=cell%ts(i)-cell%tg(i)
            del_tsl=cell%ts(i)-cell%tl(i)
            del_tgl=cell%tg(i)-cell%tl(i)
            c_ht=c_ig_hf*del_tsg+c_if_hg*del_tsl
            bm1(i3,2)= c_ht+qvol_liq(i)*dt+c_gf*del_tgl+ediff_liq(i)*dtvol &
                      -gamma_wall(i)*hi_liq_w*dt+qporous_liq(i)*dtvol
            bm1(i3,1)=-c_ht+qvol_gas(i)*dt-c_gf*del_tgl+ediff_gas(i)*dtvol &
                      +gamma_wall(i)*hi_gas_w*dt+qporous_gas(i)*dtvol      &
                      +mediff_gas(i)*dtvol
            bm1(i3,3)=-c_ig*del_tsg-c_if*del_tsl+gamma_wall(i)*dt
!
!...........Droplet continuity equation
!
            dsesd=cell%entr(i)-cell%dentr(i)
            a=-(c_ig_dtsde-c_ig_dtvde+c_if_dtsde)
            b=-(c_ig*(cell%dtsdx(i)-cell%dtgdx(i))+c_if*cell%dtsdx(i))
            c=c_ig*del_tsg+c_if*del_tsl
            sm(i3,4,1)=a*cell%yeta(i)
            sm(i3,4,2)=ad*cell%drholde(i)+c_if_dtlde*cell%yeta(i)
            sm(i3,4,3)=0.d0
            sm(i3,4,4)=cell%rhol(i)
            sm(i3,4,5)=b*cell%yeta(i)
            sm(i3,4,6)=ad*cell%drholdp(i)-c_ig*del_dtvdp*cell%yeta(i)-c_if*del_dtldp*cell%yeta(i)
            bm1(i3,4)=c*cell%yeta(i)+dsesd*dt
!
!...........Continuous liquid continuity equation
!
            yeta1=1.d0-cell%yeta(i)
            sm(i3,5,1)=a*yeta1
            sm(i3,5,2)=al*cell%drholde(i)+c_if_dtlde*yeta1
            sm(i3,5,3)=-cell%rhol(i)
            sm(i3,5,4)=-cell%rhol(i)
            sm(i3,5,5)=b*yeta1
            sm(i3,5,6)=al*cell%drholdp(i)-c_ig*del_dtvdp*yeta1-c_if*del_dtldp*yeta1
            bm1(i3,5)=c*yeta1-dsesd*dt-gamma_wall(i)*dt
! 
!...........Include the Hik gradient for an implicit treatment of HIk 
!           Gradient w.r.t temperature is taken into account
!
            cig=cig*del_tsg
            cif=cif*del_tsl
            cgf=cgf*del_tgl
            cig_hf=cig*hi_liq
            cif_hg=cif*hi_gas
            a=(cig_hf+cgf)*dHgdtg(i)*cell%dtgde(i)
            b=cif_hg*dHldtl(i)*cell%dtlde(i)
            sm(i3,1,1)=sm(i3,1,1)+a
            sm(i3,2,1)=sm(i3,2,1)-a
            sm(i3,1,2)=sm(i3,1,2)+b
            sm(i3,2,2)=sm(i3,2,2)-b
            sm(i3,3,1)=sm(i3,3,1)+cig*dHgdtg(i)*cell%dtgde(i)
            sm(i3,3,2)=sm(i3,3,2)+cif*dHldtl(i)*cell%dtlde(i)
            sm(i3,4,1)=sm(i3,4,1)-cell%yeta(i)*cig*dHgdtg(i)*cell%dtgde(i)
            sm(i3,4,2)=sm(i3,4,2)-cell%yeta(i)*cif*dHldtl(i)*cell%dtlde(i)
            sm(i3,5,1)=sm(i3,5,1)-yeta1*cig*dHgdtg(i)*cell%dtgde(i)
            sm(i3,5,2)=sm(i3,5,2)-yeta1*cif*dHldtl(i)*cell%dtlde(i)
!
!...........NCG continuity equation
!
            dtvol=dt*volpr(i)
            a=cell%alphag(i)*cell%quala(i)
            sm(i3,6,1)=a*cell%drhogde(i)
            sm(i3,6,2)=0.d0
            sm(i3,6,3)=cell%rhog(i)*cell%quala(i)
            sm(i3,6,4)=0.d0
            sm(i3,6,5)=a*cell%drhogdx(i)+MAX(1.d-10,cell%alphag(i))*cell%rhog(i)
            sm(i3,6,6)=a*cell%drhogdp(i)
            bm1(i3,6)=mdiff_gas(i)*dtvol ! Implement scalar_mass_diffusion 2015.07.29 JHLee (SNU) 
         ENDDO
!
!........Radiation Option
!
         IF(rad_model.ne.0)THEN
!!DIR$ ASSUME_ALIGNED rad_source:avx
            DO i=is,ie
               IF(npb(i).ne.0) CYCLE
               i3=i-is+1
               bm1(i3,1)=bm1(i3,1)+rad_source(i)*dtvol
            ENDDO
         ENDIF         
!
!........EVVD Option
!
         IF(l_subchannel)THEN
!DIR$ VECTOR ALIGNED
            DO i=is,ie
               IF(npb(i).ne.0) CYCLE
               i3=i-is+1
               bm1(i3,1)=bm1(i3,1)+(tm_eng_g(i)+vd_eng_g(i))*dt
               bm1(i3,2)=bm1(i3,2)+(tm_eng_l(i)+vd_eng_l(i))*dt
               bm1(i3,3)=bm1(i3,3)+(tm_mas_g(i)+vd_mas_g(i))*dt
               bm1(i3,5)=bm1(i3,5)+(tm_mas_l(i)+vd_mas_l(i))*dt
            ENDDO
         ENDIF
!
!........Mixing Vane
!
         IF(l_mixing_vane .and. time.ge.0.1d0) then
            DO i=is,ie
               IF(npb(i).ne.0) CYCLE
               i3=i-is+1
               bm1(i3,2)=bm1(i3,2)+mixing_vane_l(3,i)*dtvol
               bm1(i3,5)=bm1(i3,5)+mixing_vane_l(2,i)*dtvol
            ENDDO
         ENDIF            
!
!........Modify vapor and liquid energy equations to make thermal equilibrium for small phase fraction
!
         IF(l_th_equil.gt.0)THEN
            DO i=is,ie
               IF(npb(i).ne.0) CYCLE
               i3=i-is+1
               n_rhs=ia_nrhs(i)
               IF(cell%alphag(i).lt.alphag_min)THEN
                  sm(i3,1,:)=0.d0
                  sm(i3,1,1)=1.d0
                  bm1(i3,1)=0.d0
                  DO j=1,n_rhs
                     bm(1+j,1,i3)=0.d0
                  ENDDO
               ENDIF
               IF(cell%alphal(i).lt.alphal_min)THEN
                  sm(i3,2,:)=0.d0
                  sm(i3,2,2)=1.d0
                  bm1(i3,2)=0.d0
                  DO j=1,n_rhs
                     bm(1+  n_rhs+j,2,i3)=0.d0
                     bm(1+2*n_rhs+j,2,i3)=0.d0
                  ENDDO
               ENDIF
            ENDDO
         ENDIF
         DO i=is,ie
            IF(npb(i).ne.0) CYCLE
            i3=i-is+1
               bm(1,1,i3)=bm1(i3,1)
               bm(1,2,i3)=bm1(i3,2)
               bm(1,3,i3)=bm1(i3,3)
               bm(1,4,i3)=bm1(i3,4)
               bm(1,5,i3)=bm1(i3,5)
               bm(1,6,i3)=bm1(i3,6)
         ENDDO
!
!........Solve the 6x6 system nb blocks
!
         CALL luinverse6mn(sm,bm,nblk,1+3*ns,nb,npb(is),ia_nrhs(is))
!
!........Save the matrix elements as global variables
!
         DO i=is,ie
            IF(npb(i).ne.0) CYCLE
            i3=i-is+1
            n_rhs=ia_nrhs(i)
            DO k=1,6
               sb(i,k)=bm(1,k,i3)
            ENDDO
!
            ii=iptr_nb_k(i)
            j1=0
            IF(ii.gt.0) THEN
               j=j1
               DO ip=ia_nb(ii),ia_nb(ii+1)-1
                  j=j+1
                  sfg_non_k(ip,1)=bm(1        +j,1,i3)
                  sfl_non_k(ip,1)=bm(1+  n_rhs+j,1,i3)
                  sfd_non_k(ip,1)=bm(1+2*n_rhs+j,1,i3)
                  sfg_non_k(ip,2)=bm(1        +j,2,i3)
                  sfl_non_k(ip,2)=bm(1+  n_rhs+j,2,i3)
                  sfd_non_k(ip,2)=bm(1+2*n_rhs+j,2,i3)
                  sfg_non_k(ip,3)=bm(1        +j,3,i3)
                  sfl_non_k(ip,3)=bm(1+  n_rhs+j,3,i3)
                  sfd_non_k(ip,3)=bm(1+2*n_rhs+j,3,i3)
                  sfg_non_k(ip,4)=bm(1        +j,4,i3)
                  sfl_non_k(ip,4)=bm(1+  n_rhs+j,4,i3)
                  sfd_non_k(ip,4)=bm(1+2*n_rhs+j,4,i3)
                  sfg_non_k(ip,5)=bm(1        +j,5,i3)
                  sfl_non_k(ip,5)=bm(1+  n_rhs+j,5,i3)
                  sfd_non_k(ip,5)=bm(1+2*n_rhs+j,5,i3)
               ENDDO
               j=j1
               DO ip=ia_nb(ii),ia_nb(ii+1)-1
                  j=j+1
                  sfg6_non_k(ip)=bm(1        +j,6,i3)
                  sfl6_non_k(ip)=bm(1+  n_rhs+j,6,i3)
                  sfd6_non_k(ip)=bm(1+2*n_rhs+j,6,i3)
               ENDDO
               j1=j 
            ENDIF
            DO nf_number=0,2
               ii=ip2(nf_number)
               istart1=istart_nb1(1,nf_number)
               len    =istart_nb1(2,nf_number)
               IF(ii.gt.len) cycle
               i1=istart1+ii
!
!..............Find pointer ii to access _nb that correspond to i cell
!..............and save for next i search
!
400            CONTINUE
               IF(icell_nb(i1).lt.i) THEN
                  ii=ii+1
                  IF(ii.gt.len) goto 410
                  i1=istart1+ii
                  ip2(nf_number)=ii
                  goto 400
               ENDIF
!
               IF(icell_nb(i1).eq.i) THEN
                  istart =istart_nf(1,nf_number)
                  istart1=istart_nb1(1,nf_number)
                  i0=istart1+ii
                  j=j1
                  DO ip=ia_nb(i0),ia_nb(i0+1)-1
                     i1=istart+ip
                     j=j+1
                     sfg_nf(i1,1)=bm(1        +j,1,i3)
                     sfl_nf(i1,1)=bm(1+  n_rhs+j,1,i3)
                     sfd_nf(i1,1)=bm(1+2*n_rhs+j,1,i3)
                     sfg_nf(i1,2)=bm(1        +j,2,i3)
                     sfl_nf(i1,2)=bm(1+  n_rhs+j,2,i3)
                     sfd_nf(i1,2)=bm(1+2*n_rhs+j,2,i3)
                     sfg_nf(i1,3)=bm(1        +j,3,i3)
                     sfl_nf(i1,3)=bm(1+  n_rhs+j,3,i3)
                     sfd_nf(i1,3)=bm(1+2*n_rhs+j,3,i3)
                     sfg_nf(i1,4)=bm(1        +j,4,i3)
                     sfl_nf(i1,4)=bm(1+  n_rhs+j,4,i3)
                     sfd_nf(i1,4)=bm(1+2*n_rhs+j,4,i3)
                     sfg_nf(i1,5)=bm(1        +j,5,i3)
                     sfl_nf(i1,5)=bm(1+  n_rhs+j,5,i3)
                     sfd_nf(i1,5)=bm(1+2*n_rhs+j,5,i3)
                  ENDDO
                  j=j1
                  DO ip=ia_nb(i0),ia_nb(i0+1)-1
                     i1=istart+ip
                     j=j+1
                     sfg6_nf(i1)=bm(1        +j,6,i3)
                     sfl6_nf(i1)=bm(1+  n_rhs+j,6,i3)
                     sfd6_nf(i1)=bm(1+2*n_rhs+j,6,i3)
                  ENDDO
                  j1=j 
                  ip2(nf_number)=ip2(nf_number)+1
               ENDIF
410            CONTINUE
            ENDDO
         ENDDO
      ENDDO
!
      END SUBROUTINE scalar_matrix
