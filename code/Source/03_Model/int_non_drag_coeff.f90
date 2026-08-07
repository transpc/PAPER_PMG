!
      SUBROUTINE int_non_drag_coeff
!
!     This routine calculates coefficients for non-drag forces
!
      USE VOL_DATA                
      USE Zconst1     , ONLY: nlift,ntdf,nwlf
      USE Zflowregime , ONLY: alphag_bc,alphag_cm
      USE Zndforce    , ONLY: cwlf,clift,ctd
      USE Zvector     , ONLY: vrel_o
      USE Zzone       , ONLY: ncell_fluid
!
      IMPLICIT NONE
!
      INTEGER i
!      
      REAL(8) weight
      REAL(8) Eotvos,dh_bubble,EotvosDH,fEotvosDH,Reg
      REAL(8) gravity
!
      DATA gravity/9.81d0/
!
!
!........Tomiayma's lift coefficient model
!
       IF(nlift.eq.-2.0d0)THEN
          DO i=1,ncell_fluid
            Eotvos=gravity*(cell%rhol(i)-cell%rhog(i))*cell%D1(i)**2/cell%sigma(i)
            dh_bubble=cell%D1(i)*(1.0d0+0.163d0*Eotvos**0.757d0)**(1.d0/3.d0) ! maximum horizontal DIMENSION
            Reg=cell%rhol(i)*vrel_o(i)*cell%D1(i)/cell%lviscosl(i)
            EotvosDH=gravity*(cell%rhol(i)-cell%rhog(i))*dh_bubble**2/cell%sigma(i)
            fEotvosDH=0.00105d0*EotvosDH**3-0.0159d0*EotvosDH**2-0.0204d0*EotvosDH+0.474d0
            IF(EotvosDH.lt.4.0d0)THEN
               Clift(i)=DMIN1(0.288d0*DTANH(0.121d0*Reg),fEotvosDH)
            ELSEIF(EotvosDH.ge.4.0d0.and.EotvosDH.le.10.7d0)THEN
               Clift(i)=fEotvosDH
            ELSEIF(EotvosDH.gt.10.7d0)THEN
               Clift(i)=-0.27838d0
            ENDIF
          ENDDO
!      
!........Not use non_drag force  (-1: OFF, 0: Constant coefficient, 1: Mutplied by tuning factor)
!
       ELSEIF(nlift.eq.-1.0d0)THEN
          !Do nothing
       ELSEIF(nlift.eq.0.0d0)THEN     
          DO i=1,ncell_fluid
            Clift(i)=0.01d0
          ENDDO
       ELSE
          DO i=1,ncell_fluid
            Clift(i)=nlift 
          ENDDO
       ENDIF
!
       IF(nwlf.eq.-1.0d0)THEN
          !Do nothing
       ELSEIF(nwlf.eq.0.0d0)THEN     
          DO i=1,ncell_fluid
            Cwlf(i)=1.0d0
          ENDDO
       ELSE
          DO i=1,ncell_fluid
            Cwlf(i)=nwlf 
          ENDDO
       ENDIF      
!
       IF(ntdf.eq.-1.0d0)THEN
          !Do nothing
       ELSEIF(ntdf.eq.0.0d0)THEN     
          DO i=1,ncell_fluid
            Ctd(i)=0.1d0
          ENDDO
       ELSE
          DO i=1,ncell_fluid
            Ctd(i)=ntdf 
          ENDDO
       ENDIF
!      
!.....Relaxation for void fraction
!
      DO i=1,ncell_fluid
         IF(cell%regime(i).eq.13)THEN
            IF(nlift.ne.-1.0d0)Clift(i)=0.0d0
            IF(nwlf.ne.-1.0d0)Cwlf(i)=0.0d0   
            IF(ntdf.ne.-1.0d0)Ctd(i)=0.0d0             
         ELSEIF(cell%regime(i).eq.12)THEN
            weight=DMAX1(0.0d0,DMIN1(1.0d0,(alphag_cm-cell%alphag(i))/(alphag_cm-alphag_bc)))
            IF(nlift.ne.-1.0d0)Clift(i)=Clift(i)*weight
            IF(nwlf.ne.-1.0d0)Cwlf(i)=Cwlf(i)*weight 
            IF(ntdf.ne.-1.0d0)Ctd(i)=Ctd(i)*weight      
         ENDIF   
      ENDDO
!       
      RETURN
      END SUBROUTINE int_non_drag_coeff
