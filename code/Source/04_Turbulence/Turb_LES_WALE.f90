!
      SUBROUTINE Turb_LES_WALE(ik)
!
!     This subroutine calcultes turbulent viscosity using LES-WALE model.
!
      USE VOL_DATA
      USE Zb_condition , ONLY: vb_liq,vin_liq,vb_gas,vin_gas
      USE Zbc_index    , ONLY: icell_type
      USE Zcoord3      , ONLY: vol
      USE Zmpi         , ONLY: ncell_fp
      USE Zparam       , ONLY: ndim
      USE Zturb        , ONLY: wvis_liq,wvis_gas
      USE Zvector      , ONLY: vl_o,vg_o
      USE Zzone        , ONLY: ncell_fluid
!
      IMPLICIT NONE
!
      INTEGER ik,i
!
      REAL(8) s11,s22,s33,s12,s13,s23,sijsij
      REAL(8) g11,g22,g33,g12,g13,g21,g23,g31,g32
      REAL(8) g11_2,g22_2,g33_2,g12_2,g13_2,g21_2,g23_2,g31_2,g32_2,gkk2
      REAL(8) s11d,s22d,s33d,s12d,s13d,s23d,sijdsijd
      REAL(8) Cw,deltaV,CwdeltaV,numer,denom
      REAL(8) tviscos
      REAL(8) rho(ncell_fp),dvdx(ncell_fp,ndim,ndim)
!
      IF(ik.eq.1)THEN          !turb_phase=1 : gas only
         CALL grad_vel(1,vg_o,dvdx,vb_gas,vin_gas)
         rho(:)=cell%rhog(:)
      ELSEIF(ik.eq.2)THEN      !turb_phase=2 : liquid only
         CALL grad_vel(2,vl_o,dvdx,vb_liq,vin_liq)
         rho(:)=cell%rhol(:)
      ENDIF
!
      DO i=1,ncell_fluid
         s11=dvdx(i,1,1)
         s22=dvdx(i,2,2)
         s33=dvdx(i,3,3)
         s12=0.5d0*(dvdx(i,1,2)+dvdx(i,2,1))
         s13=0.5d0*(dvdx(i,1,3)+dvdx(i,3,1))
         s23=0.5d0*(dvdx(i,2,3)+dvdx(i,3,2))
         sijsij=s11*s11+s22*s22+s33*s33+2.0d0*(s12*s12+s13*s13+s23*s23)
!!!!!
         g11=dvdx(i,1,1)
         g22=dvdx(i,2,2)
         g33=dvdx(i,3,3)
         g12=dvdx(i,1,2)
         g13=dvdx(i,1,3)
         g21=dvdx(i,2,1)
         g23=dvdx(i,2,3)
         g31=dvdx(i,3,1)
         g32=dvdx(i,3,2)
!
         g11_2=g11*g11+g12*g21+g13*g31
         g22_2=g21*g12+g22*g22+g23*g32
         g33_2=g31*g13+g32*g23+g33*g33
         g12_2=g11*g12+g12*g22+g13*g32
         g13_2=g11*g13+g12*g23+g13*g33
         g21_2=g21*g11+g22*g21+g23*g31
         g23_2=g21*g13+g22*g23+g23*g33
         g31_2=g31*g11+g32*g21+g33*g31
         g32_2=g31*g12+g32*g22+g33*g32
         gkk2=g11_2+g22_2+g33_2
!
         s11d=g11_2-gkk2/3.0d0
         s22d=g22_2-gkk2/3.0d0
         s33d=g33_2-gkk2/3.0d0
         s12d=0.5d0*(g12_2+g21_2)
         s13d=0.5d0*(g13_2+g31_2)
         s23d=0.5d0*(g23_2+g32_2)
         sijdsijd=s11d*s11d+s22d*s22d+s33d*s33d+2.0d0*(s12d*s12d+s13d*s13d+s23d*s23d)
!!!!!
         Cw=0.5d0    !!! input
         deltaV=vol(i)**0.33333d0
         CwdeltaV=(Cw*deltaV)**2
         numer=sijdsijd**1.5d0
         denom=sijsij**2.5d0+sijdsijd**1.25d0
!
         IF(denom.eq.0.0d0)THEN
            tviscos=0.0d0
         ELSE
            tviscos=rho(i)*CwdeltaV*numer/denom
         ENDIF
!
         IF(ik.eq.1)THEN        !turb_phase=1 : gas only
            cell%tviscosg(i)=tviscos
            IF(icell_type(i).eq.1) wvis_gas(i)=cell%tviscosg(i)
            cell%eviscosg(i)=cell%lviscosg(i)+cell%tviscosg(i)
         ELSEIF(ik.eq.2)THEN    !turb_phase=2 : liquid only
            cell%tviscosl(i)=tviscos
            IF(icell_type(i).eq.1) wvis_liq(i)=cell%tviscosl(i)
            cell%eviscosl(i)=cell%lviscosl(i)+cell%tviscosl(i)
         ENDIF
      ENDDO
!
      RETURN
      END SUBROUTINE Turb_LES_WALE
!
!--------------------------------------------------------------------------------
!
      SUBROUTINE Turb_LES_cond(ik)
!
      USE VOL_DATA
      USE Zbc_index , ONLY: icell_type
      USE Zparam    , ONLY: prt
      USE Zturb     , ONLY: wvis_liq,wvis_gas,wcd_liq,wcd_gas
      USE Zzone     , ONLY: ncell_fluid
!
      IMPLICIT NONE
!
      INTEGER ik,i
!
      IF(ik.eq.1)THEN        !turb_phase=1 : gas only
         DO i=1,ncell_fluid
            IF(icell_type(i).eq.1) wcd_gas(i)=cell%lcondg(i)+wvis_gas(i)*cell%cpg(i)/prt
            cell%condg(i)=cell%lcondg(i)+cell%tviscosg(i)*cell%cpg(i)/prt
         ENDDO
      ELSEIF(ik.eq.2)THEN    !turb_phase=2 : liquid only
         DO i=1,ncell_fluid
            IF(icell_type(i).eq.1) wcd_liq(i)=cell%lcondl(i)+wvis_liq(i)*cell%cpl(i)/prt
            cell%condl(i)=cell%lcondl(i)+cell%tviscosl(i)*cell%cpl(i)/prt
         ENDDO
      ENDIF
!
      RETURN
      END SUBROUTINE Turb_LES_cond
