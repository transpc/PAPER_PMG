!
      SUBROUTINE choke_cell_avg
!
      USE VOL_DATA     , ONLY: cell
      USE Zparam       , ONLY: pi
      USE Zconst1      , ONLY: restart
      USE Zcore        , ONLY: np
      USE Zrv_choke    , ONLY:  vl_choke,vg_choke,vl_choke_o,vg_choke_o,                                 &
                                cell_leng_avg,ar_liq_avg,ar_gas_avg,p_avg,p_avg_out,pps_avg,vl_o_avg,vg_o_avg,     &
                                rhog_avg,eg_avg,el_avg,quala_avg,quals_avg,tl_avg,alphag_avg,alphal_avg, &
                                cvao_cell_avg,dcva_cell_avg,uao_cell_avg,ra_cell_avg,          &
                                icell_throat,ocell_throat,ick,num_throatface,ick_dir
      USE Zncg         , ONLY: cvao_cell,uao_cell,dcva_cell,ra_cell                                 
      USE Zvector      , ONLY: vl_o,vg_o
      USE Zcoord2      , ONLY: cell_leng    
      USE Zuserdefined , ONLY: user_rary
      USE Zpress      , ONLY: p
      USE Zzone        , only: ncell_fluid   
      USE Zmpi
      USE Ztimecon
!
      IMPLICIT NONE
!
      INTEGER :: i,ock
      LOGICAL,SAVE:: initial=.true.
!
!.....Specify the cell index with break faces
!
      cell_leng_avg=0.d0
      ar_liq_avg=0.d0
      ar_gas_avg=0.d0
      p_avg=0.d0
      p_avg_out=0.d0
      pps_avg=0.d0
      vl_o_avg=0.d0
      vg_o_avg=0.d0
      rhog_avg=0.d0
      eg_avg=0.d0
      el_avg=0.d0
      quala_avg=0.d0
      quals_avg=0.d0
      tl_avg=0.d0
      alphag_avg=0.d0
      alphal_avg=0.d0
      cvao_cell_avg=0.d0
      dcva_cell_avg=0.d0
      uao_cell_avg=0.d0
      ra_cell_avg=0.d0
!
      IF(initial) THEN
         vl_choke=0.d0
         vg_choke=0.d0
         vl_choke_o=0.d0
         vg_choke_o=0.d0   
         initial=.false.
      ENDIF
!
      DO i=1,num_throatface
         ick=icell_throat(i)
         ock=ocell_throat(i)
         
         IF(ick.gt.ncell_fluid) CYCLE
!
!.....Define the average values
!
         cell_leng_avg=cell_leng_avg+cell_leng(ick,ick_dir) 
         ar_liq_avg=ar_liq_avg+cell%alphal(ick)*cell%rhol(ick)
         ar_gas_avg=ar_gas_avg+cell%alphag(ick)*cell%rhog(ick)
         p_avg=p_avg+p(ick)                
         p_avg_out=p_avg_out+p(ock)
         pps_avg=pps_avg+cell%pps(ick)          
         vl_o_avg=vl_o_avg+vl_o(ick,ick_dir)    
         vg_o_avg=vg_o_avg+vg_o(ick,ick_dir)    
         rhog_avg=rhog_avg+cell%rhog(ick)       
         eg_avg=eg_avg+cell%eg(ick)             
         el_avg=el_avg+cell%el(ick)             
         quala_avg=quala_avg+cell%quala(ick)    
         quals_avg=quals_avg+cell%quals(ick)    
         tl_avg=tl_avg+cell%tl(ick)             
         alphag_avg=alphag_avg+cell%alphag(ick)    
         alphal_avg=alphal_avg+cell%alphal(ick)     
         cvao_cell_avg=cvao_cell_avg+cvao_cell(ick) 
         dcva_cell_avg=dcva_cell_avg+dcva_cell(ick) 
         uao_cell_avg=uao_cell_avg+uao_cell(ick)    
         ra_cell_avg=ra_cell_avg+ra_cell(ick)       
      ENDDO
!      
      IF(num_throatface.ne.0) THEN
         cell_leng_avg=cell_leng_avg/DBLE(num_throatface)
         ar_liq_avg=ar_liq_avg/DBLE(num_throatface)
         ar_gas_avg=ar_gas_avg/DBLE(num_throatface)
         p_avg=p_avg/DBLE(num_throatface)
         p_avg_out=p_avg_out/DBLE(num_throatface)
         pps_avg=pps_avg/DBLE(num_throatface)
         vl_o_avg=vl_o_avg/DBLE(num_throatface)
         vg_o_avg=vg_o_avg/DBLE(num_throatface)
         rhog_avg=rhog_avg/DBLE(num_throatface)
         eg_avg=eg_avg/DBLE(num_throatface)
         el_avg=el_avg/DBLE(num_throatface)
         quala_avg=quala_avg/DBLE(num_throatface)
         quals_avg=quals_avg/DBLE(num_throatface)
         tl_avg=tl_avg/DBLE(num_throatface)
         alphag_avg=alphag_avg/DBLE(num_throatface)
         alphal_avg=alphal_avg/DBLE(num_throatface)
         cvao_cell_avg=cvao_cell_avg/DBLE(num_throatface)
         dcva_cell_avg=dcva_cell_avg/DBLE(num_throatface)
         uao_cell_avg=uao_cell_avg/DBLE(num_throatface)
         ra_cell_avg=ra_cell_avg/DBLE(num_throatface)
      ENDIF   
!
!.....restart update       
!         
      IF(restart.ne.0)THEN
         vl_choke=user_rary(6)
         vg_choke=user_rary(7)
         vl_choke_o=user_rary(8)
         vg_choke_o=user_rary(9)            
      ENDIF   
!
      IF(np.gt.1) THEN
         CALL allreducei_r1(cell_leng_avg)  
         CALL allreducei_r1(ar_liq_avg)
         CALL allreducei_r1(ar_gas_avg)
         CALL allreducei_r1(p_avg)
         CALL allreducei_r1(p_avg_out)
         CALL allreducei_r1(pps_avg)
         CALL allreducei_r1(vl_o_avg)
         CALL allreducei_r1(vg_o_avg)
         CALL allreducei_r1(rhog_avg)
         CALL allreducei_r1(eg_avg)
         CALL allreducei_r1(el_avg)
         CALL allreducei_r1(quala_avg)
         CALL allreducei_r1(quals_avg)
         CALL allreducei_r1(tl_avg)
         CALL allreducei_r1(alphag_avg)
         CALL allreducei_r1(alphal_avg)         
         CALL allreducei_r1(cvao_cell_avg)
         CALL allreducei_r1(dcva_cell_avg)
         CALL allreducei_r1(uao_cell_avg)
         CALL allreducei_r1(ra_cell_avg)
      ENDIF   
!     
      RETURN
      END SUBROUTINE choke_cell_avg

    