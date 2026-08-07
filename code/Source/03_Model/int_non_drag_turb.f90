!
      SUBROUTINE int_non_drag_turb(Cb_bubble)
!
!.....This routine calculates turbulent dispersion force
!
      USE VOL_DATA                 
      USE Zzone        , ONLY: ncell_fluid
      USE Zconst1      , ONLY: iturb
      USE Zndforce     , ONLY: ctd
      USE Zturb        , ONLY: turb_ke
      USE Zndforce     , ONLY: s_turb_disp      
      USE Zvector      , ONLY: vrel_o
!     USE Zvoid        , ONLY: daldx
!
      IMPLICIT NONE
!
      INTEGER i
!
      REAL(8) Atdf
      REAL(8) CD,d_bubble,Pr,reciprocal_alphag,reciprocal_alphal,Coeff
      REAL(8) Cb_bubble(ncell_fluid)
!
!     not needed anymore was used previously in pressure_solve
!     CALL grad_press(cell%alphal,daldx,3)
!           IF(np.gt.1) CALL communicate_2d(daldx)
!
      DO i=1,ncell_fluid
!          
!........Use of user defined function
!    
         IF(s_turb_disp.eq.'udf')THEN
            CALL udfn_turb_disp       
!
!........CFX model
!
         ELSEIF(s_turb_disp.eq.'cfx')THEN
            IF(iturb.lt.0) STOP '### CFX-Turbulence dispersion model should be used with a proper turbulence model (iturb>=0) ###'
            ctd(i)=-Ctd(i)*cell%rhol(i)*turb_ke(i)
!
!........Gosman model
!
         ELSEIF(s_turb_disp.eq.'Gosman')THEN
            CD=Cb_bubble(i)  
            d_bubble=cell%d1(i)
            Pr=cell%eviscosl(i)/cell%lviscosl(i)
            Atdf=-3.d0/4.d0*CD*cell%tviscosl(i)/(d_bubble*Pr)
            ctd(i)=ctd(i)*Atdf*vrel_o(i)     ! delete '*cell%rhol(i)'
!
!........Burns model
!
         ELSEIF(s_turb_disp.eq.'Burns')THEN
            CD=Cb_bubble(i)
            IF(cell%alphag(i).lt.1.0d-3)THEN
               reciprocal_alphag=1.0d3
            ELSE
               reciprocal_alphag=1.0d0/cell%alphag(i)
            ENDIF
            IF(cell%alphal(i).lt.1.0d-3)THEN
               reciprocal_alphal=1.0d3
            ELSE
               reciprocal_alphal=1.0d0/cell%alphal(i)
            ENDIF
            Coeff=-Ctd(i)*CD*cell%eviscosg(i)*reciprocal_alphag/0.9d0
            !0.9=tubulent schmidt number         
            ctd(i)=Coeff*(reciprocal_alphal+reciprocal_alphag)
!			       
         ENDIF
!
      ENDDO
!
      RETURN
      END SUBROUTINE int_non_drag_turb
