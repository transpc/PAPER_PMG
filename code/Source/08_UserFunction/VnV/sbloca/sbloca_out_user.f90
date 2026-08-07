!
      SUBROUTINE sbloca_out_yhy_user
!
!     This routine writes the calculation results for 2D_loca problem
!
      USE Zinterface
      USE VOL_DATA     , ONLY: cell
      USE Zzone        , ONLY: ncell_fluid
      USE Zcore        , ONLY: np
      USE Zvec_param   , ONLY: nf_non,nf_mcc,nf_inl,nf_out,nf_flux
      USE Znum_cell    , ONLY: i_neigh,indexr_sort,                      &
                               nf_number_nb,lens,nf_number_id,istart_nfs
      USE Zconst2      , ONLY: dt
      USE Ztimecon     , ONLY: time
      USE Zb_condition , ONLY: eb_liq,rhob_liq
      USE Zcoord3      , ONLY: vol
      USE Zpress       , ONLY: p
      USE Zqvol        , ONLY: h_il,h_ig,gamma,qvol_liq
      USE Zsbloca      , ONLY: break_flow,break_flow_eng_l,break_flow_eng_g,si_flow,si_flow_eng, &
                               q_liq,break_flow_eng,break_flow_int,break_flow_eng_int,           &
                               break_flow_eng_l_int,break_flow_eng_g_int,si_flow_int,            &
                               si_flow_eng_int,q_liq_int,ncell_si,ge_err,eng_gg,pw_l,            &
                               pw_g,ge_err_int,eng_gg_int,pw_g_int
      USE Zvec_major   , ONLY: flux_l_nf,flux_g_nf,flux_d_nf, &
                               void_conv_nf
!
      IMPLICIT NONE
!      
!.....Local variables
      INTEGER :: i,k,k0,k1
      INTEGER :: i1
      REAL(8) :: hi_gas,hi_liq,t_s,t_l,t_g,err,gg,avg
      REAL(8) :: flux_l,flux_d,flux_g
!.....Local arrays
      REAL(8) :: tmp(6)
      REAL(8) :: aavg(ncell_fluid)
!.....Local vector arrays
      REAL(8) :: avg_nf(nf_flux)
!      
      break_flow      =0.d0
      break_flow_eng_g=0.d0
      break_flow_eng_l=0.d0
      si_flow         =0.d0
      si_flow_eng     =0.d0
!
      IF(time.gt.300.d0)THEN
         IF(ncell_si(1).gt.0)THEN
!...........carefull (2,k) must fall within non,inl,out,mcc ONLY
!...........neighbor have been sorted use indexr_sort to retrieve j=2   
            k=ncell_si(1)
            k0=i_neigh(k)-1
            k1=indexr_sort(2+k0)
!
!...........Get offset i1 in vector space of (k1,k)
!
            CALL get_vector_disp(k1,k,i1)
            IF(i1.gt.0) THEN
               flux_l=flux_l_nf(i1)
               flux_d=flux_d_nf(i1)
               flux_g=flux_g_nf(i1)
            ELSEIF(i1.lt.0) THEN
               i1=-i1
               flux_l=-flux_l_nf(i1)
               flux_d=-flux_d_nf(i1)
               flux_g=-flux_g_nf(i1)
            ENDIF
           
!
            break_flow      = cell%alphal(k)*cell%rhol(k)           *flux_l &
                             +cell%alphad(k)*cell%rhol(k)           *flux_d &
                             +cell%alphag(k)*cell%rhog(k)           *flux_g
            break_flow_eng_g= cell%alphag(k)*cell%rhog(k)*cell%eg(k)*flux_g
            break_flow_eng_l= cell%alphal(k)*cell%rhol(k)*cell%el(k)*flux_l &
                             +cell%alphad(k)*cell%rhol(k)*cell%el(k)*flux_d
         ENDIF
!
         IF(ncell_si(2).gt.0)THEN
!...........carefull (2,k) must fall within non,inl,out,mcc ONLY
!...........neighbor have been sorted use indexr_sort to retrieve j=2   
            k=ncell_si(2)
            k0=i_neigh(k)-1
            k1=indexr_sort(2+k0)
!
!...........Get offset i1 in vector space of (k1,k)
!
            CALL get_vector_disp(k1,k,i1)
            IF(i1.gt.0) THEN
               flux_l=flux_l_nf(i1)
            ELSEIF(i1.lt.0) THEN
               i1=-i1
               flux_l=-flux_l_nf(i1)
            ENDIF
!
            si_flow    =-rhob_liq(2)          *flux_l
            si_flow_eng=-rhob_liq(2)*eb_liq(2)*flux_l
         ENDIF
!
         q_liq=0.d0
         DO i=1,ncell_fluid
            q_liq=q_liq+qvol_liq(i)*vol(i)
         ENDDO
!
         IF(np.gt.1)THEN
            tmp(1)=q_liq
            tmp(2)=break_flow
            tmp(3)=break_flow_eng_g
            tmp(4)=break_flow_eng_l
            tmp(5)=si_flow
            tmp(6)=si_flow_eng
            CALL allreducei_r(tmp,6)
            q_liq           =tmp(1)
            break_flow      =tmp(2)
            break_flow_eng_g=tmp(3)
            break_flow_eng_l=tmp(4)
            si_flow         =tmp(5)
            si_flow_eng     =tmp(6)
         ENDIF
!
         break_flow_eng      =break_flow_eng_g    +break_flow_eng_l
         break_flow_int      =break_flow_int      +break_flow      *dt
         break_flow_eng_int  =break_flow_eng_int  +break_flow_eng  *dt
         break_flow_eng_g_int=break_flow_eng_g_int+break_flow_eng_g*dt
         break_flow_eng_l_int=break_flow_eng_l_int+break_flow_eng_l*dt
         si_flow_int    =si_flow_int    +si_flow*dt
         si_flow_eng_int=si_flow_eng_int+si_flow_eng*dt
         q_liq_int      =q_liq_int      +q_liq      *dt
!
      ENDIF
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
      DO i1=1,nf_flux
         avg_nf(i1)=void_conv_nf(i1)*flux_g_nf(i1)
      ENDDO
      CALL sum_nf(0,-1,        &
                  avg_nf,aavg)
!
      ge_err=0.d0
      eng_gg=0.d0
      pw_g=0.d0
      pw_l=0.d0
      DO i=1,ncell_fluid
         IF(gamma(i).ge.0.d0)THEN
            hi_gas=cell%hgsat(i)
            hi_liq=cell%hl(i)
         ELSE
            hi_gas=cell%hg(i)
            hi_liq=cell%hlsat(i)
         ENDIF
         t_s=cell%ts(i)
         t_g=cell%tg(i)
         t_l=cell%tl(i)
         err=( H_il(i)*(t_s-t_l)+H_ig(i)*(t_s-t_g)  &
              +gamma(i)*(hi_gas-hi_liq))*vol(i)*dt
         ge_err=ge_err+err
         gg=-vol(i)*( hi_liq*H_ig(i)*(cell%ts(i)-cell%tg(i))         &
                     +hi_gas*H_il(i)*(cell%ts(i)-cell%tl(i)))/(hi_gas-hi_liq)
         eng_gg=eng_gg+gg*dt
!
         avg =aavg(i)
         pw_g=pw_g-p(i)*avg*dt
!
      ENDDO
!
      IF(np.gt.1)THEN
         tmp(1)=ge_err
         tmp(2)=eng_gg
         tmp(3)=pw_g
         CALL allreducei_r(tmp,3)
         ge_err=tmp(1)
         eng_gg=tmp(2)
         pw_g  =tmp(3)
      ENDIF
!
      ge_err_int=ge_err_int+ge_err*dt
      eng_gg_int=eng_gg_int+eng_gg
      pw_g_int  =pw_g_int  +pw_g
!
      END SUBROUTINE sbloca_out_yhy_user
!----------------------------------------------------------------------
!
      SUBROUTINE sbloca_out_user
!
!     Save output for TECPLOT
!
      USE VOL_DATA , ONLY: cell
      USE Zzone    , ONLY: ncell_fluid                            
      USE Zcore    , ONLY: np,myrank      
      USE Ztimecon , ONLY: time,cfl_ratio
      USE Zare     , ONLY: are_liq,are_gas,are_drp
      USE Zcoord3  , ONLY: vol
      USE Zsbloca  , ONLY: p_break,t_break,q_break,h_break,a_break,si_flow,si_flow_int,si_flow_eng_int, &
                            rhom_break,break_flow_int,break_flow_eng_int,break_flow_eng_l_int,           &
                            break_flow_eng_g_int,q_liq_int,q_liq,eng_gg,eng_gg_int,ge_err,ge_err_int,    &
                            pw_l,pw_g,pw_g_int,break_flow,break_flow_eng,si_flow_eng
!
      IMPLICIT NONE 
!
!.....Local variables
      INTEGER i
      LOGICAL,SAVE::  INITIAL      
      REAL(8),SAVE :: rcs_mass0, rcs_energy0
      REAL(8) rcs_liq_vol, rcs_mass,rcs_energy,rcs_mass_s,rcs_energy_s,err_mass, err_eng
      REAL(8) rcs_energy_g, rcs_energy_l,rcs_energy_s_g,rcs_energy_s_l
!.....Local arrays
      REAL(8) :: tmp(11)
!
      DATA INITIAL /.TRUE./      
!
!.....SBLOCA output
!
      rcs_liq_vol =0.d0
      rcs_mass    =0.d0
      rcs_energy  =0.d0
      rcs_energy_g=0.d0
      rcs_energy_l=0.d0
      DO i=1,ncell_fluid
         rcs_liq_vol=rcs_liq_vol+vol(i)*cell%alphal(i)
         rcs_mass=rcs_mass+vol(i)*cell%rhom(i)
         rcs_energy_g=rcs_energy_g+vol(i)*are_gas(i)
         rcs_energy_l=rcs_energy_l+vol(i)*(are_liq(i)+are_drp(i))
         rcs_energy=rcs_energy_g+rcs_energy_l 
      END DO
!
      IF(np.gt.1)THEN
         tmp( 1)=rcs_liq_vol
         tmp( 2)=rcs_mass
         tmp( 3)=rcs_energy_g
         tmp( 4)=rcs_energy_l
         tmp( 5)=rcs_energy
         tmp( 6)=p_break
         tmp( 7)=t_break
         tmp( 8)=q_break
         tmp( 9)=rhom_break
         tmp(10)=h_break
         tmp(11)=a_break
         CALL allreducei_r(tmp,11)
         rcs_liq_vol =tmp( 1)
         rcs_mass    =tmp( 2)
         rcs_energy_g=tmp( 3)
         rcs_energy_l=tmp( 4)
         rcs_energy  =tmp( 5)
         p_break     =tmp( 6)
         t_break     =tmp( 7)
         q_break     =tmp( 8)
         rhom_break  =tmp( 9)
         h_break     =tmp(10)
         a_break     =tmp(11)
      ENDIF
!
      IF(INITIAL)THEN
         rcs_mass0=rcs_mass
         rcs_energy0=rcs_energy
         INITIAL=.FALSE.
      ENDIF
      rcs_mass_s = rcs_mass0 - break_flow_int + si_flow_int
      rcs_energy_s = rcs_energy0 - break_flow_eng_int + si_flow_eng_int + q_liq_int + pw_g_int
      rcs_energy_s_g = - break_flow_eng_g_int + eng_gg_int + pw_g_int
      rcs_energy_s_l = rcs_energy0 - break_flow_eng_l_int + si_flow_eng_int + q_liq_int - eng_gg_int
      err_mass = (rcs_mass_s-rcs_mass) / rcs_mass
      err_eng = (rcs_energy_s-rcs_energy) / rcs_energy
!
      IF(myrank.eq.0)THEN
         WRITE(20,100) time, cfl_ratio, p_break, t_break, q_break, rhom_break, h_break, a_break,   &
                       break_flow, break_flow_int, si_flow, si_flow_int, rcs_liq_vol, rcs_mass,    &
                       rcs_mass_s, q_liq, q_liq_int, ge_err, ge_err_int, break_flow_eng,           &
                       break_flow_eng_int, si_flow_eng, si_flow_eng_int, rcs_energy, rcs_energy_s, &
                       err_mass, err_eng,rcs_energy_g,rcs_energy_s_g,rcs_energy_l,rcs_energy_s_l,  &
                       eng_gg,eng_gg_int,pw_g,pw_l
      ENDIF
100   FORMAT(35(e15.8,1x))
!
      END SUBROUTINE sbloca_out_user
