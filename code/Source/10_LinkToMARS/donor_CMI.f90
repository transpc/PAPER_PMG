!
      SUBROUTINE donor_CMI(time_inst)
!
!     Set the donor properties at the CUPID-MARS interface.
!        c3dpv(i,1): alphag*rhog*eg
!        c3dpv(i,2): alphal*rhol*el
!        c3dpv(i,3): alphag*rhog
!        c3dpv(i,4): alphad*rhol
!        c3dpv(i,5): alphal*rhol
!        c3dpv(i,6): alphag*rhog*x
!        c3dpv(i,7): alphad*rhol*el!                                                                      
!     c3rtp(i,j,k) contains the interfacing 1D cell data.
!        k=1: liquid-phase volume fraction,
!        k=2: vapor-phase volume fraction (including noncondensible gases
!        k=3: noncondensible gases fraction (X),
!        k=4: liquid-phase density,
!        k=5: vapor-phase density,
!        k=6: liquid-phase internal energy,
!        k=7: vapor-phase internal energy.
!        k=8: entrained droplet/vapor interfacial area density,
!        k=9: pressure.
!
      USE vol_data
      USE Zparam
      USE c3com_cupid, only:c3dpv,mcdirect,mcgdirect,i3cupid
      USE Zcore        , ONLY: cupid_alone,cupid_mars,myrank,np      
      USE Zmars        , ONLY: ncupvol,ncell_old
      USE Zb_condition , ONLY: vb_liq,vb_drp,vb_gas,vin_liq,vin_drp,vin_gas,&
                                  alphab_liq,alphab_gas,alphab_drp, &
                                  rhob_liq, rhob_gas, rhob_drp,&
                                  eb_liq, eb_gas, eb_drp
      USE Zzone        , ONLY: ncell_fluid
      USE c3com_cupid  , ONLY: i3invtbl
      USE Zbc_index    , ONLY: vin_norm      
      USE Ztimecon     , ONLY: time
!
      IMPLICIT NONE
!
      INCLUDE 'c3com.h'
!
      LOGICAL,SAVE:: initial     
!      
      INTEGER(4) i, ii
      INTEGER(4) i_mb,j,idx,k
      INTEGER(4),SAVE,ALLOCATABLE :: nbc_cup(:),nbc_cup_ngh(:),icell_cup(:)
!
      REAL(8),SAVE,ALLOCATABLE :: vliq_cup(:),vliq_cup_o(:)
      REAL(8),SAVE :: changetime,deltime
!      
      real(8) time_inst,a2_gas 
! 
      DATA initial/.true./      
! 
!DEC$IF defined (SPACE)
      CALL c3com_copy_S2C
!DEC$ENDIF  
!      
      CALL broadcast_r(c3rtp(1,:,1),72)
      CALL broadcast_r(c3rtp(1,:,2),72)
      CALL broadcast_r(c3rtp(1,:,3),72)
      CALL broadcast_r(c3rtp(1,:,4),72)
      CALL broadcast_r(c3rtp(1,:,5),72)
      CALL broadcast_r(c3rtp(1,:,6),72)
      CALL broadcast_r(c3rtp(1,:,7),72)
      CALL broadcast_r(c3rtp(1,:,10),72)
!      
      CALL broadcast_r(c3vl(1,:),72)
      CALL broadcast_r(c3vg(1,:),72)
!      
      CALL broadcast_r(c3pa(1,:),72)     !p
      CALL broadcast_r(c3delp(1,:),72)   !dp(pressure correction)
      CALL broadcast_r(c3alphf(1,:),72)  !vl*
      CALL broadcast_r(c3alphg(1,:),72)  !vg*
      CALL broadcast_r(c3betaf(1,:),72)  !betal-geometrical pressure coefficient
      CALL broadcast_r(c3betag(1,:),72)  !betag-
      CALL broadcast_r(c3area(1,:),72)   !convection area
!
      CALL broadcast_r(c3xi(1,:),72)       
      DO i=1,72 
         CALL broadcast_r(c3yeta(1,i,:),72)
      ENDDO
!         
      DO i=1,i3nic(2)
!          
         c3dpv(i,1)=c3rtp(1,i,2)*c3rtp(1,i,5)*c3rtp(1,i,7) !alphag*rhog*eg
         c3dpv(i,2)=c3rtp(1,i,1)*c3rtp(1,i,4)*c3rtp(1,i,6) !alphal*rhol*el
         c3dpv(i,3)=c3rtp(1,i,2)*c3rtp(1,i,5)              !alphag*rhog
         c3dpv(i,4)=0.0d0                                  !alphad*rhol      
         c3dpv(i,5)=c3rtp(1,i,1)*c3rtp(1,i,4)              !alphal*rhol
         c3dpv(i,6)=c3rtp(1,i,2)*c3rtp(1,i,5)*c3rtp(1,i,3) !alphag*rhog*x
         c3dpv(i,7)=0.0d0                                  !alphad*rhol*el                                 
         c3dpv(i,8)=c3rtp(1,i,1)*c3rtp(1,i,4)*c3rtp(1,i,10) !alphal*rhol*cboron                 
!
!......direction for donor property of convective terms
!                                                                       cupid mars
         IF(c3vl(1,i).lt.0.0d0)mcdirect(i)=-1 !from mars to cupid  -|v| in    out
         IF(c3vl(1,i).ge.0.0d0)mcdirect(i)=1  !from cupid to mars  +|v| out   in
         IF(c3vg(1,i).lt.0.0d0)mcgdirect(i)=-1
         IF(c3vg(1,i).ge.0.0d0)mcgdirect(i)=1
!         mcdirect(i)=i3dir(1,i)   direction of volume connection
!         mcgdirect(i)=i3dir(1,i)  -1=from mars out to cupid in, 1=from cupid out to mars in
!         
      ENDDO
!
      RETURN
      END SUBROUTINE donor_CMI
