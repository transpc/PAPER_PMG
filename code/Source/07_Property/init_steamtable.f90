!
      SUBROUTINE init_steamtable
!
!     This routine initializes various thermodynamic properties at the boundaries
!
      USE STM_TBL_cupid , ONLY: st_tbl,          &
                                nt,ndxstd,nfluid
      USE Zbc_index    , ONLY: nvin,npin
      USE Zb_condition , ONLY: alphab_gas,alphab_liq,alphab_drp,        &
                               alpha_gas_nd,alpha_liq_nd,alpha_drp_nd,  &
                               e_gas_nd,e_liq_nd,e_drp_nd,              &
                               eb_gas,eb_liq,eb_drp,                    &
                               rho_gas_nd,rho_liq_nd,rho_drp_nd,        &
                               rhob_gas,rhob_liq,rhob_drp,              &
                               t_gas_nd,t_liq_nd,t_drp_nd,              &
                               tb_gas,tb_liq,tb_drp,                    &
                               qualab,quala_nd,pbnd,p_fb,               &
                               lvisb_liq,lvisb_gas
      USE Zncg         , ONLY: tao,cvao_nvin,uao_nvin,dcva_nvin,ra_nvin,             &
                                cvao_npin,uao_npin,dcva_npin,ra_npin,qn_nvin,qn_npin
!
      IMPLICIT NONE 
!
!.....Local variables
      INTEGER :: i
      LOGICAL :: erx
      REAL(8) :: tmp1,tmp2
!.....Local arrays
      REAL(8) :: s(36),ts_b(1)
      REAL(8) :: qn_npin0(8),qn_nvin0(8)
!
!.....Flow boundary properties at cell face
!
      DO i=1, nvin
         alphab_liq(i) = 1.d0 - alphab_gas(i)
         alphab_drp(i) = 0.d0
      ENDDO
!      
      DO i=1,nvin
!      
!........liquid, steam/gas      
!
         qn_nvin0(:)=qn_nvin(i,:) 
         CALL convert_temp2erg(p_fb(i),tb_liq(i),tb_gas(i),qualab(i),eb_liq(i),eb_gas(i),rhob_liq(i),rhob_gas(i),tmp1,tmp2, &
                               tao,cvao_nvin(i),uao_nvin(i),dcva_nvin(i),ra_nvin(i),qn_nvin0)
!
!........Droplet property is same as liquid
!
         tb_drp(i)   = tb_liq(i)
         rhob_drp(i) = rhob_liq(i)
         eb_drp(i)   = eb_liq(i)         
!
!........viscosity
!
         s(2)=p_fb(i)
         s(9)=0.d0
         IF(nfluid.eq.1)then 
            CALL sth2x2_cupid(st_tbl(ndxstd),st_tbl(ndxstd+nt),s,erx)
         ELSEIF(nfluid.eq.2)then 
            CALL std2x2_cupid(st_tbl(ndxstd),s,erx)
         ELSEIF(nfluid.eq.15)then 
            CALL nth2x2_cupid(st_tbl(ndxstd),s,erx) 
         ELSE 
            CALL strpx_cupid(st_tbl(ndxstd),s,erx) 
         ENDIF           
         ts_b(1)=s(1)
         CALL viscos_lw_cupid(tb_gas(i),rhob_gas(i),ts_b,lvisb_gas(i),'vap',1)
         CALL viscos_lw_cupid(tb_liq(i),rhob_liq(i),ts_b,lvisb_liq(i),'liq',1)
      ENDDO
!
!.....Pressure boundary properties at cell center
!
      DO i=1, npin
         alpha_liq_nd(i) = 1.d0 - alpha_gas_nd(i)
         alpha_drp_nd(i) = alpha_liq_nd(i)
      ENDDO
!
      DO i=1,npin
!      
!........liquid, steam/gas      
!      
         qn_npin0(:)=qn_npin(i,:) 
         CALL convert_temp2erg(pbnd(i),t_liq_nd(i),t_gas_nd(i),quala_nd(i),e_liq_nd(i),       &
                               e_gas_nd(i),rho_liq_nd(i),rho_gas_nd(i),tmp1,tmp2,             &
                               tao,cvao_npin(i),uao_npin(i),dcva_npin(i),ra_npin(i),qn_npin0)
!
!........Droplet property is same as liquid
!
         t_drp_nd(i)   = t_liq_nd(i)
         rho_drp_nd(i) = rho_liq_nd(i)
         e_drp_nd(i)   = e_liq_nd(i)         
      ENDDO
!
      END SUBROUTINE init_steamtable

