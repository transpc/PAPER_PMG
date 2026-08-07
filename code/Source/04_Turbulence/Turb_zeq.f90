!
      SUBROUTINE Turb_zeq
!
!     This routine calculates turbulent viscosities base upon
!     zero equation model.
!
      USE VOL_DATA                 
      USE Zzone        , ONLY: ncell_fluid
      USE Zparam       , ONLY: ndim
      USE Zcoord1      , ONLY: xloc
      USE Zmpi         , ONLY: ncell_fp
      USE Zndforce     , ONLY: d_bfc,dvdxl
      USE Zturb        , ONLY: utau,yplus
      USE Zturbzeq     , ONLY: cell_hindex,chheight,tleng,vorticity
      USE Zturbzeq     , ONLY: s_turb_zero
      USE Zvector      , ONLY: vg_n,vl_n
!
      IMPLICIT NONE
!
      INTEGER i,j
      INTEGER nChHeight   
!      
      LOGICAL,SAVE :: INITIAL,INITIAL1      
!      
      REAL(8) cmu_zeq,te_frac,k_en
      REAL(8) fmu,Ut
      REAL(8) Dee
      REAL(8) wvorticity
      REAL(8) Y,Ymax, LMIXVORTmax,Fmax,UdIF,Fwake1,Fwake2,Fwake,Fkleb
      REAL(8),PARAMETER::kappa=0.41d0,alpha=0.0168d0,Ckleb=0.3d0,Cwk=1.0d0,Ccp=1.6d0
!
      DATA INITIAL/.TRUE./
      DATA INITIAL1 /.TRUE./
!      
      REAL(8), SAVE, ALLOCATABLE:: Lmixvort(:),vlabs(:)
!      
      IF(initial1) THEN 
            ALLOCATE(Lmixvort(ncell_fp),vlabs(ncell_fp))
            initial1=.false.
      ENDIF      
!
      DO i=1,ncell_fluid
         yplus(i)=d_bfc(i)*utau(i)*cell%rhol(i)/cell%lviscosl(i)
      ENDDO
!
      SELECT CASE(s_turb_zero)
!
!.....Constant model
!      
      CASE('constant')
         DO i=1,ncell_fluid
            cell%tviscosl(i)=cell%lviscosl(i) * 100.d0
            cell%tviscosg(i)=cell%lviscosg(i) * 100.d0
         ENDDO
!
!.....Telluride model
!
      CASE('telluride')
         DO i=1,ncell_fluid
            cmu_zeq=0.05d0
            te_frac=0.01d0
            k_en=te_frac * 0.5D0*dot_product(vl_n(i,:),vl_n(i,:))
            cell%tviscosl(i)=cell%rhol(i)*cmu_zeq*dsqrt(k_en)*tleng(i)
            k_en=te_frac * 0.5D0*dot_product(vg_n(i,:),vg_n(i,:))
            cell%tviscosg(i)=cell%rhog(i)*cmu_zeq*dsqrt(k_en)*tleng(i)
         ENDDO
!
!.....CFX model
!
      CASE('cfx')
         DO i=1,ncell_fluid
            fmu=0.005d0
            Ut=dot_product(vl_n(i,:),vl_n(i,:))
            cell%tviscosl(i)=cell%rhol(i)*fmu*dsqrt(Ut)*tleng(i)
            Ut=dot_product(vg_n(i,:),vg_n(i,:))
            cell%tviscosg(i)=cell%rhog(i)*fmu*dsqrt(Ut)*tleng(i)
         ENDDO
!
!.....Prandtl model
!         
      CASE('prandtl')
         DO i=1,ncell_fluid
            wvorticity=(dvdxl(i,1,2)-dvdxl(i,2,1))**2
            If(ndim.eq.3)wvorticity=wvorticity+(dvdxl(i,1,ndim)-dvdxl(i,ndim,1))**2 &
                                              +(dvdxl(i,2,ndim)-dvdxl(i,ndim,2))**2
            vorticity(i)=wvorticity**0.5d0
         ENDDO
!
         DO i=1,ncell_fluid
             tleng(i)=kappa*d_bfc(i)
             cell%tviscosl(i)=cell%rhol(i)*tleng(i)*tleng(i)*vorticity(i)
         ENDDO         
!
!.....Constant model
!
      CASE('vanDriest')
          DO i=1,ncell_fluid
              wvorticity=(dvdxl(i,1,2)-dvdxl(i,2,1))**2
              If(ndim.eq.3)wvorticity=wvorticity+(dvdxl(i,1,ndim)-dvdxl(i,ndim,1))**2 &
                                                +(dvdxl(i,2,ndim)-dvdxl(i,ndim,2))**2
              vorticity(i)=wvorticity**0.5d0
          ENDDO
!
         DO i=1,ncell_fluid
            tleng(i)=kappa*d_bfc(i)*Dee(yplus(i))
            cell%tviscosl(i)=cell%rhol(i)*tleng(i)*tleng(i)*vorticity(i)
         ENDDO
!
!...Baldwin_Lomax
!
      CASE('baldwin_Lomax')
         IF(initial)THEN
            initial=.false.
            nChHeight=100
            ChHeight(1,0)=MINVAL(xloc(:,ndim))
            ChHeight(1,nChHeight)=MAXVAL(xloc(:,ndim))
            DO i=1,nChHeight-1
               ChHeight(1,i)=(ChHeight(1,nChHeight)-ChHeight(1,0))/nChHeight*i !channel ChHeight from bottom
            ENDDO
            DO i=1,ncell_fluid
               DO j=1,nChHeight
                  IF(xloc(i,ndim).gt.ChHeight(1,j-1) .and. xloc(i,ndim).le.ChHeight(1,j))&
                          cell_Hindex(i)=j !remember  index of icell's ChHeight
               ENDDO
            ENDDO
         ENDIF
!
         DO i=1,ncell_fluid
            wvorticity=(dvdxl(i,1,2)-dvdxl(i,2,1))**2
            If(ndim.eq.3)wvorticity=wvorticity+(dvdxl(i,1,ndim)-dvdxl(i,ndim,1))**2 &
                                              +(dvdxl(i,2,ndim)-dvdxl(i,ndim,2))**2
            wvorticity=wvorticity**0.5d0
            vorticity(i)=wvorticity
            tleng(i)=kappa*d_bfc(i)*Dee(yplus(i)) !¸Â³ª?
            LmixVort(i)=tleng(i)*wvorticity
            vlabs(i)=dsqrt(dot_product(vl_n(i,:),vl_n(i,:)))
         ENDDO
!
         ChHeight(2:5,:)=0.d0
         DO i=1,ncell_fluid
            IF(ChHeight(2,cell_Hindex(i)).lt.LmixVort(i))THEN 
               ChHeight(2,cell_Hindex(i))=LmixVort(i) 
               ChHeight(3,cell_Hindex(i))=d_bfc(i)   
               ChHeight(4,cell_Hindex(i))=dsqrt(dot_product(vl_n(i,:),vl_n(i,:)))
            ENDIF
            IF(ChHeight(5,cell_Hindex(i)).lt.vlabs(i))THEN
               ChHeight(5,cell_Hindex(i))=vlabs(i)
            ENDIF
         ENDDO
         ChHeight(4,:)=ChHeight(5,:)-ChHeight(4,:) 
!         
         DO i=1,ncell_fluid
            IF(yplus(i).le.600)THEN 
               tleng(i)=kappa*d_bfc(i)*Dee(yplus(i))
               cell%tviscosl(i)=cell%rhol(i)*tleng(i)*tleng(i)*vorticity(i)
            ELSE
               Y=d_bfc(i) 
               Ymax=ChHeight(3,cell_Hindex(i))
               LMIXVORTmax=ChHeight(2,cell_Hindex(i))
               Fmax=1.0d0/0.41d0*LMIXVORTmax
               UdIF=ChHeight(4,cell_Hindex(i))
               Fwake1=Ymax*Fmax
               Fwake2=Fwake1
               IF(Fmax.ne.0.d0)Fwake2=1.0*Ymax*UdIF**2./Fmax
               Fwake=dmin1(Fwake1,Fwake2)
               Fkleb=1.0d0/(1.0d0+(Y/(Ymax/Ckleb))**6.)
               cell%tviscosl(i)=alpha*Ccp*Fwake*Fkleb
            ENDIF              
         ENDDO
!
!.....Noto model
!        
      CASE('noto')
         DO i=1,ncell_fluid
            cell%tviscosl(i)=0.4d0*cell%lviscosl(i)*yplus(i)*(1.0d0-exp(-0.0017d0*(yplus(i)**2)))
            cell%tviscosl(i)=dmin1(5000.0d0*cell%lviscosl(i),cell%tviscosl(i))
         ENDDO
!
!...User model
!        
      CASE('user')
         CALL udfn_turb_zeq
      END SELECT
!
      DO i=1,ncell_fluid
         tleng(i)=tleng(i)*utau(i)*cell%rhol(i)/cell%lviscosl(i)
      ENDDO
!
      RETURN
      END SUBROUTINE Turb_zeq
!
      Function Dee(yplus)
!     
      IMPLICIT NONE
!      
      REAL(8) yplus, Dee
      REAL(8), PARAMETER:: aplus=26.d0
!      
      Dee=1.0d0-dexp(-yplus/aplus)
 !    
      END FUNCTION Dee
