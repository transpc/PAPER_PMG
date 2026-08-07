!-------------------------------------------------------------------------------------- -
COMMON / c3com / i_where, i3nic, i3cell, &
flag_relap, flag_cobra, &
dt_super, dt_relap, dt_cobra, &
c3pa, c3uf, c3ug, c3al, c3arxq, c3rhof, c3rhog, c3brn, c3vpgno, &
c3rtp, c3alphg, c3alphf, c3betag, c3betaf, c3xi, c3yeta, c3area, &
c3vg, c3vl, c3delp

LOGICAL flag_relap     !whether MARS calculaton is succeeded or not
LOGICAL flag_cobra     !whether CUPID calculaton is succeeded or not
REAL(8) dt_super       !dominating time step
REAL(8) dt_relap       !MARS time step
REAL(8) dt_cobra       !CUPUD time step
!CUPVOLs(CUPID->MARS) properties(MARS->CUPID) at interface
REAL(8) c3pa(1, 72)     !pressure in CUPVOL, calculated in CUPID
REAL(8) c3ug(1, 72)     !internal energy of gas in CUPVOL, calculated in CUPID
REAL(8) c3uf(1, 72)     !internal energy of liquid in CUPVOL, calculated in CUPID
REAL(8) c3al(1, 72)     !gas volume fraction in CUPVOL, calculated in CUPID
REAL(8) c3arxq(1, 72)   !NCG mass fraction in CUPVOL, calculated in CUPID
REAL(8) c3rhof(1, 72)   !density of liquid in CUPVOL, calculated in CUPID
REAL(8) c3rhog(1, 72)   !density of gas in CUPVOL, calculated in CUPID
REAL(8) c3brn(1, 72)    !boron concentration in CUPVOL, calculated in CUPID
REAL(8) c3vpgno(1, 72)  !vapor generation in CUPVOL, currently zero, calculated in CUPID
REAL(8) c3rtp(1, 72, 15) !convective properties in MARS 1D cell, which is neighbor of CUPVOL(see donor_CMI.f90), calculated in MARS
!c3rtp(1, i, 1) = donor liquid volume fraction
!c3rtp(1, i, 2) = donor gas volume fraction
!c3rtp(1, i, 3) = donor gas NCG mass fraction
!c3rtp(1, i, 4) = donor liquid density
!c3rtp(1, i, 5) = donor gas density
!c3rtp(1, i, 6) = donor liquid internal energy
!c3rtp(1, i, 7) = donor gas internal energy
!c3rtp(1, i, 10) = donor boron concentration !mass ratio
!pressure coefficients, dPr(1:n) = c3xi(1:n) + c3yeta(1, 1:n, 1 : n)dPc(1:n), Vc = c3alphf + c3betaf * d(dP)
REAL(8) c3alphg(1, 72)   !temporal velocity of gas, calculated in MARS
REAL(8) c3alphf(1, 72)   !temporal velocity of liquid, calculated in MARS
REAL(8) c3betag(1, 72)   !geometrical pressure gradient coefficient of gas, calculated in MARS
REAL(8) c3betaf(1, 72)   !geometrical pressure gradient coefficient of gas, calculated in MARS
REAL(8) c3xi(1, 72)      !pressure correction of interface MARS cell, calculated in MARS
REAL(8) c3yeta(1, 72, 72) !mutiplier to pressure correction of CUPVOL cells, calculated in MARS
REAL(8) c3area(1, 72)    !convective area between CUPVOL and MARS 1D cell, read in MARS
!new velocity, pressureand pressure correction
REAL(8) c3vg(1, 72)     !updated velocity at interface, calculated in CUPID
REAL(8) c3vl(1, 72)     !updated velocity at interface, calculated in CUPID
REAL(8) c3delp(1, 72)   !pressure correction in CUPVOL, calculated in CUPID
INTEGER i_where        !calculation step identifier
INTEGER i3nic(3)       !i3nic(2) means the number of CUPVOLSand i3nic(1) = 0, i3nic(3) = i3nic(2) + i3nic(1)
INTEGER i3cell(1, 72)
!----------------------------------------------------------------------------------
COMMON / c3com_dll / c3time_sys, c3tend, & !r1
c3rktpow_ctl_val, c3kfactor_ctl_val,c3RPV3d_ctl_val, & !r2
c3mflow_junleg, & !r3
c3p_volleg, c3tl_volleg,c3p_tmdpvol2nd,c3t_tmdpvol2nd,c3run_mode3_dur,c3run_mode3_dpcri,c3run_mode2_dur,c3run_mode2_dpcri, & !r4
i3myrank, i3np, i3cupid_alone, i3cupid_mars, &
i3cplmaster, i3cplmars, i3marsin, &
i3rx_trip, i3rcp_trip, i3degba_trip, i3mslb_trip, i3rod_trip, &
i3n_junleg, i3n_volleg, i3run_mode, i3n_tmdpvol2nd
REAL(8) c3time_sys   !time, defined in SYSTEM
REAL(8) c3tend       !end time, defined in CUPID
INTEGER i3myrank     !core index, defined in SYSTEM(0~np - 1)
INTEGER i3np         !the number of cores, defined in SYSTEM
INTEGER i3cupid_alone!indicator of CUPID alone, defined in SYSTEM(0, 1)
INTEGER i3cupid_mars !indicator of coupling CUPID / SYSTEM, defined in SYSTEM(0, 1)
INTEGER i3cplmaster  !stage of coupling CUPID / NEUTRONICS, defined in CUPID(0, 1)
INTEGER i3cplmars    !stage of coupling CUPID / SYSTEM, defined in CUPID(0, 1, 2, 3)
INTEGER i3marsin     !
INTEGER i3rx_trip    !status of reactor trip, defined in SYSTEM(0, 1)
INTEGER i3rcp_trip   !status of RCP trip, defined in SYSTEM(0, 1)
INTEGER i3degba_trip !status of DEGBA trip, defined in SYSTEM(0, 1)
INTEGER i3mslb_trip  !status of MSLB trip, defined in SYSTEM(0, 1)
INTEGER i3rod_trip   !status of control rod trip, defined in SYSTEM(0, 1)
REAL(8) c3rktpow_ctl_val  !power ratio to normal power
REAL(8) c3kfactor_ctl_val !kfactor of 1D RCP volume to match pressure loss of 1D and 3D RPV
REAL(8) c3RPV3d_ctl_val   !3D RPV pressure drop
INTEGER i3n_junleg        !the number of junctions
INTEGER i3n_volleg        !the number of volumes
REAL(8) c3mflow_junleg(20)!mass flows at the junctions
REAL(8) c3p_volleg(20)    !liquid temperatures at the volumes
REAL(8) c3tl_volleg(20)	  !pressures at the volumes
REAL(8) c3p_tmdpvol2nd(2) !pressure at tmdpvol at 2nd stage
REAL(8) c3t_tmdpvol2nd(2) !temperature at tmpdvol at 2nd stage
INTEGER i3run_mode        !if i3run_mode == 2, generate indta_1st_correction.i; 3, indta_2nd_correction.i
REAL(8) c3run_mode3_dur,c3run_mode3_dpcri,c3run_mode2_dur,c3run_mode2_dpcri
      INTEGER i3n_tmdpvol2nd    !the number of tmdpvol at 2nd step
!-----------------------------------------------------------------------
     COMMON/c3com_old/&
         c3pa_o,c3ug_o,c3uf_o,c3al_o,c3arxq_o,c3rhof_o,c3rhog_o,c3brn_o,c3vpgno_o,&       
         c3rtp_o,& 
         c3alphg_o,c3alphf_o,c3betag_o,c3betaf_o,c3xi_o,c3yeta_o,c3area_o,&        
         c3vg_o,c3vl_o,c3delp_o      
      REAL(8)     c3pa_o(1,72) !pressure in CUPVOL ,calculated in CUPID                                                                    
      REAL(8)     c3ug_o(1,72) !internal energy of gas in CUPVOL ,calculated in CUPID                                                      
      REAL(8)     c3uf_o(1,72) !internal energy of liquid in CUPVOL ,calculated in CUPID                                                   
      REAL(8)     c3al_o(1,72) !gas volume fraction in CUPVOL ,calculated in CUPID                                                         
      REAL(8)   c3arxq_o(1,72) !NCG mass fraction in CUPVOL ,calculated in CUPID                                                           
      REAL(8)   c3rhof_o(1,72) !density of liquid in CUPVOL ,calculated in CUPID                                                           
      REAL(8)   c3rhog_o(1,72) !density of gas in CUPVOL ,calculated in CUPID                                                              
      REAL(8)    c3brn_o(1,72) !boron concentration in CUPVOL ,calculated in CUPID                                                         
      REAL(8)  c3vpgno_o(1,72) !vapor generation in CUPVOL, currently zero ,calculated in CUPID                                            
      REAL(8) c3rtp_o(1,72,15) !convective properties in MARS 1D cell, which is neighbor of CUPVOL(see donor_CMI.f90) ,calculated in MARS   
      REAL(8)  c3alphg_o(1,72) !temporal velocity of gas, calculated in MARS                         
      REAL(8)  c3alphf_o(1,72) !temporal velocity of liquid, calculated in MARS                      
      REAL(8)  c3betag_o(1,72) !geometrical pressure gradient coefficient of gas, calculated in MARS 
      REAL(8)  c3betaf_o(1,72) !geometrical pressure gradient coefficient of gas, calculated in MARS 
      REAL(8)     c3xi_o(1,72) !pressure correction of interface MARS cell, calculated in MARS       
      REAL(8)c3yeta_o(1,72,72) !mutiplier to pressure correction of CUPVOL cells, calculated in MARS 
      REAL(8)   c3area_o(1,72) !interface area, calculated in MARS                 
      REAL(8)     c3vg_o(1,72) !updated velocity at interface ,calculated in CUPID
      REAL(8)     c3vl_o(1,72) !updated velocity at interface ,calculated in CUPID
      REAL(8)   c3delp_o(1,72) !pressure correction in CUPVOL ,calculated in CUPID    	
!----------------------------------------------------------------------------------
      COMMON/c3com_mar/ n_repet,nstep_c,nstep_r,&
         i3dir,i1Cvoln,i1Cvodn,i1Cvndx,i1jndx,i1Rvodn,i1Rvndx,& !MARS
         flag_stop,unit3_i,unit3_o,flag_bd,flag_wp,flag_bd2!not handling
      !not handling      
      LOGICAL flag_stop     !whether MARS calculaton will be stop or not
      LOGICAL unit3_i       !indicate si unit if true. use only si unit
      LOGICAL unit3_o       !indicate si unit if true. use only si unit
      LOGICAL flag_bd       !CUPVOLs exist if true. defined, but not used
      LOGICAL flag_bd2      !CUPVOLs exist if true. defined, but not used
      LOGICAL flag_wp       !CUPVOLS exist if true. defined, but not used
      !MARS
      INTEGER n_repet       !the number of calculation at this step
      INTEGER nstep_c       !the number of CUPID calculations, not used
      INTEGER nstep_r       !the number of MARS, not used
      INTEGER i3dir(1,72)   !flow direction at interface, -1:from MARS to CUPID, 1:from CUPID to MARS
      INTEGER i1Cvoln(1,72) !volume number of CUPVOLS
      INTEGER i1Cvodn(1,72) !volume index of CUPVOLS, i1Cvodn(1,k)       
      INTEGER i1jndx(1,72)  !junction index of MARS junction at interface
      INTEGER i1Rvndx(1,72) !volume index of MARS volume at interface   
      INTEGER i1Rvodn(1,72) !sourcp index at MARS system matrix
      INTEGER i1Cvndx(1,72) !not used
!----------------------------------------------------------------------------------
      COMMON/c3com_cobra/i3nodr,&!COBRA/MARS
         i3chan,i3mode,i3modet,iend50,i3line,vpackflag,hpackflag,i3bcn,ibeg50, &!COBRA
         sdbvolpack,sdbvolpacko,c3pack, &!COBRA&MARS
         overcorrection, & !COBRA
	     packfactor, & !COBRA/MARS
         c3dpmt,c3vvl,c3vvg,c3odr,c3dpm,r3ent,r3liq,cc3rtp,aloold !COBRA
      
      LOGICAL sdbvolpack
      LOGICAL sdbvolpacko(1,72)
      LOGICAL overcorrection(1,72)
      LOGICAL c3pack(100,100)  
      
      REAL(8) packfactor     
      REAL(8) cc3rtp(1,72,15)                                      
      REAL(8) r3ent(1,72)
      REAL(8) r3liq(1,72)                
      REAL(8) c3dpmt(1,72,6)
      REAL(8) c3vvl(1,72)
      REAL(8) c3vvg(1,72)
      REAL(8) c3odr(6)
      REAL(8) c3dpm(6)    
!     REAL(8) packdvf(1,72)
      REAL(8) aloold(100,100)
!
      INTEGER i3nodr
      INTEGER ibeg50
      INTEGER iend50
      INTEGER i3chan(1,72)
      INTEGER i3mode(1,72)
      INTEGER i3line(30)  
      INTEGER i3modet(1,72,3)
      INTEGER i3bcn(1,72)
      INTEGER hpackflag(100,100)
      INTEGER vpackflag(100,100)        	
