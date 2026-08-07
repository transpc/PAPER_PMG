!
      SUBROUTINE udfn_wall_lub_i(i,d_wall,Awlf)
!
!     This routine includes user defined wall lubrication force model
!
      USE VOL_DATA                 
      USE Zndforce        , ONLY: cwlf
!
      IMPLICIT NONE
!
!.....Input
      REAL(8) d_wall
!.....Output
      REAL(8) Awlf
!.....Local variables
      INTEGER i
      REAL(8) dpipe,gravity,rholi,rhogi,d_bubblei,sigmai,lviscosli,eotvos,ctomi
      REAL(8) d_bubble,vrela2,eotvosdh,dh_bubble
!   
      Dpipe=0.15d0 ! 2D: 0.2, 3D : 0.15 Marine:0.1, hydraulic-diameter !7.839514d-3 for PSBT
      gravity=9.81d0
      rholi=cell%rhol(i)
      rhogi=cell%rhog(i)
      d_bubblei=cell%D1(i)
      sigmai=cell%sigma(i)
      lviscosli=cell%lviscosl(i)
      Eotvos=gravity*(rholi-rhogi)*d_bubblei**2/sigmai
      dh_bubble=d_bubblei*(1.0d0+0.163d0*Eotvos**0.757d0)**(1.d0/3.d0) ! maximum horizontal dimension
      EotvosDH=gravity*(rholi-rhogi)*dh_bubble**2/sigmai
      IF(EotvosDH.lt.1.0d0)then
         Ctomi=0.47d0
      ELSEIF(EotvosDH.ge.1.0d0.and.EotvosDH.le.5.0d0)then
         Ctomi=EXP(-0.933d0*EotvosDH+0.179)
      ELSEIF(EotvosDH.gt.5.0d0.and.EotvosDH.le.33.0d0)then
         Ctomi=0.00599*EotvosDH-0.0187
      ELSEIF(EotvosDH.gt.33.0d0)then
         Ctomi=0.179d0
      ENDIF
      Awlf= cell%alphag(i)*cell%rhol(i)*Cwlf(i)*Vrela2*Ctomi &
           *MAX(0.0d0,d_bubble/2.0d0*(1/d_wall**2.0d0-1/(Dpipe-d_wall)**2.0d0))
!
      END SUBROUTINE udfn_wall_lub_i
!
      SUBROUTINE udfn_wall_lub_v(d_wall,Awlf)
!
!     This routine includes user defined wall lubrication force model
!
      USE VOL_DATA                 
      USE Zzone           , ONLY: ncell_fluid
      USE Zndforce        , ONLY: cwlf
!
      IMPLICIT NONE
!
!.....Input
      REAL(8) d_wall(ncell_fluid)
!.....Output
      REAL(8) Awlf(ncell_fluid)
!.....Local variables
      INTEGER i
      REAL(8) dpipe,gravity,rholi,rhogi,d_bubblei,sigmai,lviscosli,eotvos,ctomi
      REAL(8) d_bubble,vrela2,eotvosdh,dh_bubble
!   
      Dpipe=0.15d0 ! 2D: 0.2, 3D : 0.15 Marine:0.1, hydraulic-diameter !7.839514d-3 for PSBT
      gravity=9.81d0
      DO i=1,ncell_fluid
         rholi=cell%rhol(i)
         rhogi=cell%rhog(i)
         d_bubblei=cell%D1(i)
         sigmai=cell%sigma(i)
         lviscosli=cell%lviscosl(i)
         Eotvos=gravity*(rholi-rhogi)*d_bubblei**2/sigmai
         dh_bubble=d_bubblei*(1.0d0+0.163d0*Eotvos**0.757d0)**(1.d0/3.d0) ! maximum horizontal dimension
         EotvosDH=gravity*(rholi-rhogi)*dh_bubble**2/sigmai
         IF(EotvosDH.lt.1.0d0)then
            Ctomi=0.47d0
         ELSEIF(EotvosDH.ge.1.0d0.and.EotvosDH.le.5.0d0)then
            Ctomi=EXP(-0.933d0*EotvosDH+0.179)
         ELSEIF(EotvosDH.gt.5.0d0.and.EotvosDH.le.33.0d0)then
            Ctomi=0.00599*EotvosDH-0.0187
         ELSEIF(EotvosDH.gt.33.0d0)then
            Ctomi=0.179d0
         ENDIF
         Awlf(i)= cell%alphag(i)*cell%rhol(i)*Cwlf(i)*Vrela2*Ctomi &
                 *MAX(0.0d0,d_bubble/2.0d0*(1/d_wall(i)**2-1/(Dpipe-d_wall(i))**2))
      ENDDO
!
      END SUBROUTINE udfn_wall_lub_v
