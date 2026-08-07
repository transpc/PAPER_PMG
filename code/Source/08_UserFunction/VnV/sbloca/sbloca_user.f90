!
      SUBROUTINE sbloca_user
!
!      This routine declares flow properties at pressure outlet boundaries
!
!
      USE VOL_DATA      , ONLY: cell   
      USE Zzone         , ONLY: ncell_fluid 
      USE Znum_cell     , ONLY: i_neigh,indexr_sort
      USE Zbc_index     , ONLY: nbcon,nvin,vin_norm
      USE STM_TBL_cupid , ONLY: st_tbl,             &
                                nt,np,ns,ns2,ndxstd
      USE Zb_condition  , ONLY: cb_pl,cb_pg,alphab_liq,alphab_gas,rhob_liq,rhob_gas, &
                                eb_liq,eb_gas,tb_liq,tb_gas,vb_liq,vb_gas,p_fb
      USE Zcoord1       , ONLY: xloc
      USE Zpress        , ONLY: p
      USE Zpress_coeff  , ONLY: coefp_l,coefp_g
      USE Zsbloca       , ONLY: ncell_si,p_break,t_break,q_break,rhom_break,h_break, &
                                a_break,rl_break,rg_break,el_break,eg_break,cl_break,cg_break
      USE Zvec_geo      , ONLY: sa_nf
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i,k,i1,k0,k1,j0
      INTEGER :: ntab_si,ipos,ichok_mb
      INTEGER :: it    
      LOGICAL, SAVE::INITIAL       
      LOGICAL :: erx            
      REAL(8) :: tsat,psats,vsubfs,vsubgs,usubfs,usubgs,hsubfs,hsubgs, &
                 betafs,betags,kapafs,kapags,cpfs,cpgs,entfs,entgs
      REAL(8),SAVE:: flux_break, flux_si, break_area, hg_break, vg_break, flux_break_o     
!.....Local arrays
      REAL(8) s(36)
      REAL(8) xtab_si(22), ytab_si(22)
!
      DATA ntab_si / 22 /
      DATA xtab_si / 0.0d0, 0.5d6, 1.0d6, 1.5d6, 2.0d6, 2.5d6, 3.0d6, 3.5d6, 4.0d6, 4.5d6,   &
                     5.0d6, 5.5d6, 6.0d6, 6.5d6, 7.0d6, 7.5d6, 8.0d6, 8.5d6, 9.0d6, 9.5d6,   &
                     10.0d6, 10.1d6 /
      DATA ytab_si / 70.2d0, 68.6d0, 66.8d0, 65.2d0, 63.4d0, 61.5d0, 59.5d0, 57.5d0, 55.4d0, 53.4d0, &
                     51.1d0, 48.6d0, 45.9d0, 43.1d0, 40.2d0, 37.0d0, 33.4d0, 29.3d0, 24.3d0, 18.2d0, &
                     7.9d0,  0.0d0 /
      DATA INITIAL /.TRUE./
!     
      EQUIVALENCE(s( 1),tsat),    &
                 (s(10),psats),   &
                 (s(11),vsubfs),  &
                 (s(12),vsubgs),  &
                 (s(13),usubfs),  &
                 (s(14),usubgs),  &
                 (s(15),hsubfs),  &
                 (s(16),hsubgs),  &
                 (s(17),betafs),  &
                 (s(18),betags),  &
                 (s(19),kapafs),  &
                 (s(20),kapags),  &
                 (s(21),cpfs),    &
                 (s(22),cpgs),    &
                 (s(25),entfs),   &
                 (s(26),entgs)
!
!.....Turn off the pressure boundary
!
      IF(INITIAL) CALL nbcon_change_start
!
      IF(INITIAL)THEN
         ncell_si(1)=0
         ncell_si(2)=0
         DO i=1,ncell_fluid
            IF(xloc(i,1).eq.0.9d0.and.xloc(i,2).eq.9.9d0) ncell_si(1)=i
            IF(xloc(i,1).eq.2.9d0.and.xloc(i,2).eq.9.9d0) ncell_si(2)=i
         ENDDO
         flux_break=0.d0
         ichok_mb=0
      ENDIF
!
!.....Define Break
!
      p_break=0.0d0
      t_break=0.0d0
      q_break=0.0d0
      rhom_break=0.0d0
      h_break=0.0d0
      a_break=0.0d0
!
      IF(ncell_si(1).gt.0)THEN
!
         break_area=2.026829916d-3                                     ! 2 inch break
!
!........Break node properties
!
         k=ncell_si(1)
!
         p_break=p(k)
         t_break=cell%quals(k)*cell%tg(k)+(1.0d0-cell%quals(k))*cell%tl(k)
         q_break=cell%quals(k)
         rhom_break=cell%rhom(k)
         h_break=cell%quals(k)*cell%hg(k)+(1.0d0-cell%quals(k))*cell%hl(k)
         a_break=cell%alphag(k)
         rg_break=cell%rhog(k)
         rl_break=cell%rhol(k)
         eg_break=cell%eg(k)
         el_break=cell%el(k)
         cg_break=coefp_g(k)
         cl_break=coefp_l(k)
         hg_break=cell%hg(k)
         vg_break=1.d0/cell%rhog(k)
!
!........Murdock_Bauman Model
!         
         flux_break_o=flux_break
         CALL CritFlow_for_LOCA(p_break,1.d5,2.0d0/rhom_break,h_break,1,0,flux_break,ichok_mb)
!        IF (ichok_mb.eq.0) flux_break=0.1d0*flux_break+0.9d0*flux_break_o
!
!........Set break b.c.
!
         IF (ichok_mb.gt.0)THEN
            k0=i_neigh(k)-1
            k1=indexr_sort(2+k0)
!
!...........Get offset i1 in vector space of (k1,k)
!
            CALL get_vector_disp(k1,k,i1)
            i1=abs(i1)
!
            flux_break=break_area/sa_nf(i1)*flux_break
!            npb(k)=0
            nvin=1
            vin_norm(nvin)=0
            j0=i_neigh(k)-1
            nbcon(k1+j0)=nvin
            cb_pl(nvin)=coefp_l(k)
            cb_pg(nvin)=coefp_g(k)
            alphab_liq(nvin)=cell%alphal(k)
            alphab_gas(nvin)=cell%alphag(k)
            rhob_liq(nvin)=cell%rhol(k)
            rhob_gas(nvin)=cell%rhog(k)
            eb_liq(nvin)=cell%el(k)
            eb_gas(nvin)=cell%eg(k)
            tb_liq(nvin)=cell%tl(k)
            tb_gas(nvin)=cell%tg(k)
            vb_liq(nvin,2)=flux_break/cell%rhom(k)
            vb_gas(nvin,2)=vb_liq(nvin,2)
         ELSE
!            npb(k) = 1
!            nbcon(2,k) = nin_max+1
!            pbnd(1) = p(k)
!            cell(k)%alphag=1.d0
!            cell(k)%alphag_o=1.d0
!            cell(k)%alphal=0.d0
!            cell(k)%alphal_o=0.d0
         ENDIF
      ENDIF         
!
!.....Define SI
!
      IF(ncell_si(2).gt.0)THEN
         k=ncell_si(2)
!
         ipos=1
         CALL INTERP1(ntab_si,p(k),xtab_si,ytab_si,ipos,flux_si)
!
!........Set SI b.c.
!
         nvin=2
         k0=i_neigh(k)-1
         k1=indexr_sort(2+k0)
!
!........Get offset i1 in vector space of (k1,k)
!
         CALL get_vector_disp(k1,k,i1)
         i1=abs(i1)
!
         vb_liq(nvin,2)=-flux_si/sa_nf(i1)/3600.d0
         vb_gas(nvin,2)=vb_liq(nvin,2)
         j0=i_neigh(k)-1
!        nbcon(2,k)=nvin
         nbcon(k1+j0)=nvin
         cb_pl(nvin)=coefp_l(k)
         cb_pg(nvin)=coefp_g(k)
         alphab_liq(nvin)=1.0d0
         alphab_gas(nvin)=0.0d0
         tb_liq(nvin)=373.15d0
         tb_gas(nvin)=373.15d0
         p_fb(nvin)=p(k)
!
!........liquid
!
!  Initialize s for sth2x3_cupid
!        s(:)=0.d0      
         s(1)=tb_liq(nvin)
         s(2)=p_fb(nvin)
         CALL sth2x3_cupid(s,it,erx,                          &
                           st_tbl(ndxstd),                    &
                           st_tbl(ndxstd+nt),                 &
                           st_tbl(ndxstd+nt+np+13*ns+13*ns2))
         IF(erx)then
            print *, '#### ERROR: sth2x3_cupid called from sbloca_user'
            pause
            stop
         ENDIF
         rhob_liq(nvin)=1.d0/vsubfs
         eb_liq(nvin)=usubfs
!
!........steam
!
         s(2)=p_fb(nvin)
         CALL sth2x2_cupid(st_tbl(ndxstd),st_tbl(ndxstd+nt),s,erx)
         rhob_gas(nvin)=1.0d0/vsubgs
         eb_gas(nvin)=usubgs
!
!         rhob_liq(nvin) = 962.8247d0        ! 10 MPa, 100 C
!         rhob_gas(nvin) = 0.59d0
!         eb_liq(nvin) = 420000.0d0
!         eb_gas(nvin) = 2676009.0d0
!
      ENDIF
!      
      IF(INITIAL) CALL nbcon_change_end
         INITIAL=.FALSE.
!
      RETURN
      END SUBROUTINE sbloca_user
!
!------------------------------------------------------------------------------------------------
!
      SUBROUTINE INTERP1(N,X,XTBL,YTBL,IPOS0,Y)
!
!     Single table linear interpolation routine.
!     No extrapolation: If X is out of range, the table's END value will be
!
      IMPLICIT NONE
!
      INTEGER::IPOS   
      INTEGER, INTENT(IN)::N
!
!.....INITIAL guess for location (IF 0, XTBL is equi-spac, If non-zero, returns the actual location) 
!      
      INTEGER, INTENT(INOUT)::IPOS0
!
      REAL, INTENT(IN)::X
      REAL, INTENT(OUT)::Y
      REAL, INTENT(IN), DIMENSION(N)::XTBL,YTBL
      REAL :: FACT
! 
      IF(N>0)THEN
         CALL LOOKUP(N,X,XTBL,IPOS0,IPOS)
         IF(IPOS>0)THEN
            FACT=(X-XTBL(IPOS))/(XTBL(IPOS+1)-XTBL(IPOS))
            Y=YTBL(IPOS)+FACT*(YTBL(IPOS+1)-YTBL(IPOS) )
         ELSE
!         
!...........X is out of range.
!
            Y=YTBL(IABS(IPOS))
         ENDIF
      ENDIF
!      
      RETURN
      END SUBROUTINE INTERP1
!
!------------------------------------------------------------------------------------------------
!
      SUBROUTINE LOOKUP(N,X,XTBL,IPOS0,LOC)
!      
      IMPLICIT NONE
!
!     Table look-up.  Locates value X in table XTBL.
!
      INTEGER, INTENT(IN)::N
      INTEGER, INTENT(INOUT)::IPOS0              ! INITIAL guess for location (IF 0, table is equi-spac,If non-zero, returns the actual location)
      INTEGER, INTENT(OUT)::LOC                  ! Located position, Returns negative IF X is out of range
!      
      REAL, INTENT(IN)::X
      REAL, INTENT(IN), DIMENSION(N)::XTBL                                    
!
!      LOC=IMIN0(IPOS0,N-1)
      LOC=MIN0(IPOS0,N-1)
      IF((XTBL(N)-X)*(X-XTBL(1))<0.d0)THEN
!
!........Check for X being out of range, to prevent extrapolation.
!
         LOC=-N
         IF((XTBL(2)-XTBL(1))*(XTBL(1)-X)>0.0) LOC=-1
         IF(IPOS0>0) IPOS0=-LOC
!
      ELSEIF(IPOS0<=0)THEN
!
!........Calculate position in equi-spaced XTBL:
!
         LOC=IFIX((X-XTBL(1))/(XTBL(2)-XTBL(1)))+1
         LOC=MIN(N-1,LOC)
!
      ELSEIF(XTBL(2)>XTBL(1))THEN
!
!........Find position in vari-spaced XTBL, monotoniCALLy increasing.
!
         IF(X<XTBL(LOC))THEN
            LOC=LOC-1
            DO WHILE(X<XTBL(LOC))
               LOC=LOC-1
            ENDDO
         ELSE
            DO WHILE(X>XTBL(LOC+1))
               LOC=LOC+1
            ENDDO
         ENDIF
         IPOS0=LOC
!
      ELSE
!     
!........Find position in vari-spaced XTBL, monotoniCALLy decreasing.
!
         IF(X>XTBL(LOC))THEN
            LOC=LOC-1
            DO WHILE(X>XTBL(LOC))
               LOC=LOC-1
            ENDDO
         ELSE
            DO WHILE(X<XTBL(LOC+1))
               LOC=LOC+1
            ENDDO
         ENDIF
         IPOS0=LOC
      ENDIF
!      
      END SUBROUTINE LOOKUP 
!
!------------------------------------------------------------------------------------------------
!
      SUBROUTINE CritFlow_for_LOCA(press_up_in,press_down_in,NODE_SVOL,NODE_ENTH,LKFLAG,model,FLUX,ichok)
!
!     CRITNONC calculates the critical flow limit for any flow path.
!
      IMPLICIT NONE
!   
      INTEGER, INTENT(IN)  :: LKFLAG         ! Cue for leak path [0/1=momentum/leak, 2=Bernouli, 3=W(P)]
                                               ! If LKFLAG .EQ. 2, FLUX is also subcooling (input).
      INTEGER, INTENT(IN)  :: model          ! Flag for crit flow model: 0 = HEM,  1 = H-F,  2 = Moody, 3 = LOCA
      INTEGER, INTENT(OUT) :: ichok          ! Flag for choking: 0 = under chiking, 1 = choking
!            
      REAL,INTENT(OUT) :: flux           ! Critical mass flux or leak mass flux (output).   
      REAL,INTENT(IN)  :: press_up_in    ! Upstream 압력
      REAL,INTENT(IN)  :: press_down_in  ! Downstream 압력
      REAL,INTENT(IN)  :: NODE_SVOL   
      REAL press_up, press_down
      REAL NODE_ENTH
      REAL p_throat
      REAL, PARAMETER :: P_CRIT = 22.064E6
!
      IF(press_up_in<=press_down_in)THEN
         flux=0.0d0
         ichok=1
!
      ELSE
!
!........Limit of tables
!
         IF(press_up_in>=P_CRIT)THEN
            press_up=P_CRIT
            press_down=press_up-(press_up_in-press_down_in)
         ELSE
            press_up=press_up_in
            press_down=press_down_in
         ENDIF
         flux=0.0d0
         ichok=1
         CALL Get_Water_Crit_Flow_Data_from_Table(press_up,NODE_ENTH,MODEL,flux,p_throat)
         IF(LKFLAG==1.AND.p_throat<press_down)THEN
            ichok=1
            flux=DSQRT(2.0d0*(press_up-press_down)/NODE_SVOL)
         ENDIF
!      
      ENDIF
!      
      END SUBROUTINE CritFlow_for_LOCA
!
!------------------------------------------------------------------------------------------------
!
      SUBROUTINE Murdock_Bauman(p_tot,p_down,v_gas,cp_stm,LKFLAG,flux,p_throat)
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN)  :: LKFLAG       ! Cue for leak path [0/1=momentum/leak, 2=Bernouli, 3=W(P)]
                                             ! If LKFLAG .EQ. 2, FLUX is also subcooling (input).
      REAL, INTENT(IN) :: p_tot, p_down, v_gas, cp_stm
      REAL, INTENT(OUT) :: flux, p_throat
      REAL :: cv_stm
      REAL :: gamma, aa, bb, cc, dd, ee, p_ratio
!
      cv_stm=cp_stm-461.526d0
      gamma=cp_stm
      gamma=gamma/cv_stm
      aa=2.0d0/(gamma+1.0d0)
      bb=gamma/(gamma-1.0d0)
      p_throat=aa**bb*p_tot
!
      IF(p_throat<=p_down.AND.LKFLAG==1)THEN
         p_ratio=p_down/p_tot
         dd=p_ratio**(2.0d0/gamma)
         ee=1.0d0-p_ratio**(1.0d0-1.0d0/gamma)
         flux=DSQRT(2.0d0*bb*dd*ee*p_tot/v_gas)
      ELSE
         cc=(gamma+1.0d0)/(2.0d0*(gamma-1.0d0))
         flux=aa**cc*DSQRT(gamma*p_tot/v_gas)
      ENDIF
!
      END SUBROUTINE Murdock_Bauman
!
!------------------------------------------------------------------------------------------------
!
      SUBROUTINE Get_Water_Crit_Flow_Data_from_Table(p_in,h_in,model_option,mflux,p_throat)
!      
      USE Critical_Flow_Table
!   
      IMPLICIT NONE
!   
      INTEGER, INTENT(IN) :: model_option
      INTEGER :: p_index = 1
      INTEGER :: h_index = 1
!
      REAL, INTENT(IN) :: p_in,h_in   
      REAL, INTENT(OUT) :: mflux,p_throat   
!
      SELECT CASE( model_option )
      CASE (1)
         CALL TWOINTV(ENTH_DIM, PRESS_DIM, p_in, h_in, p_table, h_table, mf_hnf_table, p_index, h_index, mflux)
         CALL TWOINTV(ENTH_DIM, PRESS_DIM, p_in, h_in, p_table, h_table, pc_hnf_table, p_index, h_index, p_throat)
      CASE (2)
         CALL TWOINTV(ENTH_DIM, PRESS_DIM, p_in, h_in, p_table, h_table, mf_mdy_table, p_index, h_index, mflux)
         CALL TWOINTV(ENTH_DIM, PRESS_DIM, p_in, h_in, p_table, h_table, pc_mdy_table, p_index, h_index, p_throat)
      CASE (3)
         CALL TWOINTV(ENTH_DIM, PRESS_DIM, p_in, h_in, p_table, h_table, mf_loc_table, p_index, h_index, mflux)
         CALL TWOINTV(ENTH_DIM, PRESS_DIM, p_in, h_in, p_table, h_table, pc_loc_table, p_index, h_index, p_throat)
      CASE DEFAULT
         CALL TWOINTV(ENTH_DIM, PRESS_DIM, p_in, h_in, p_table, h_table, mf_hem_table, p_index, h_index, mflux)
         CALL TWOINTV(ENTH_DIM, PRESS_DIM, p_in, h_in, p_table, h_table, pc_hem_table, p_index, h_index, p_throat)
      ENDSELECT

      END SUBROUTINE Get_Water_Crit_Flow_Data_from_Table
!
!------------------------------------------------------------------------------------------------
!
      SUBROUTINE TWOINTV(NY,NX,X,Y,XTBL,YTBL,ZTBL,XPOS0,YPOS,Z)
!
!     Table look-up for an output variable as a function of 2 variables with linear interpolation.
!     No extrapolation: If X is out of range, the table's END value will be USEd.      
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: NX, NY
      INTEGER, INTENT(INOUT) :: XPOS0, YPOS   
      INTEGER :: POS, POSP1
!
      REAL, INTENT(IN) :: X,Y
      REAL, INTENT(IN), DIMENSION(NX) :: XTBL
      REAL, INTENT(IN), DIMENSION(NY, NX) :: YTBL,ZTBL
      REAL, INTENT(OUT) :: Z
      REAL :: FACT1,FACT2,Y1,Y2
!
      CALL LOOKUP(NX,X,XTBL,XPOS0,POS)
!
      IF (POS.GT.0.d0)THEN
         POSP1=POS+1
         FACT2=(X-XTBL(POS))/(XTBL(POSP1)-XTBL(POS))
      ELSE
!
!.... X is out of range.
!
         POS=IABS(POS)
         POSP1=POS
         FACT2=0.d0
      ENDIF
!   
      FACT1=1.d0-FACT2
!      YPOS=IMIN0(IMAX0(YPOS,1),NY-1)
      YPOS=MIN0(MAX0(YPOS,1),NY-1)
      Y1=YB(1)
      Y2=YB(NY)
      IF((Y2-Y)*(Y-Y1).LT.0.d0)THEN
!
!.....Y is out of range.
!
         IF((Y2-Y1)*(Y1-Y).GT.0.d0) YPOS=1
         Z=ZB(YPOS)
!      
         RETURN
      ENDIF
!
     IF(Y2.GT.Y1)THEN
        IF(Y.LT.YB(YPOS))THEN
           YPOS=YPOS-1
           DO WHILE(Y<YB(YPOS))
              YPOS=YPOS-1
           ENDDO
        ELSE 
           DO WHILE(Y>YB(YPOS+1))
              YPOS=YPOS+1
           ENDDO
        ENDIF
!
     ELSE
        IF(Y.GT.YB(YPOS))THEN
           YPOS=YPOS-1
           DO WHILE(Y>YB(YPOS))
              YPOS=YPOS-1
           ENDDO
        ELSE
           DO WHILE(Y<YB(YPOS+1))
              YPOS=YPOS+1
           ENDDO
        ENDIF
     ENDIF
!
     Y1=YB(YPOS)
     Y2=YB(YPOS+1)
     Z=((Y2-Y)*ZB(YPOS)+(Y-Y1)*ZB(YPOS+1))/(Y2-Y1)
     RETURN
!
     CONTAINS
!
        REAL FUNCTION YB (J)
           IMPLICIT NONE
           INTEGER, INTENT(IN) :: J
           YB=FACT1*YTBL(J,POS)+FACT2*YTBL(J,POSP1)
        END FUNCTION YB
!
        REAL FUNCTION ZB (J)
           IMPLICIT NONE
           INTEGER, INTENT(IN) :: J
           ZB=FACT1*ZTBL(J,POS)+FACT2*ZTBL(J,POSP1)
        END FUNCTION ZB
!
     END SUBROUTINE TWOINTV
