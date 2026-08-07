!
      SUBROUTINE int_area
!
!     This routine assign topology and calculates interfacial area according 
!     to the assigned topology.
!     11=bubbly,12=transition,3=sharp    gamma_void  3  3   3
!     21=transition,22=transition,3=sharp            21 22  23  
!     31=mist,32=transition,3=sharp                  11 12  13 
!                                                 0.0  bc cm  g: void fraction  
!     See SUBROUTINE initialize_topology(itp) for how to assign topology.
!
      USE VOL_DATA     , ONLY: cell
      USE Zzone        , ONLY: ncell_fluid
      USE Zconst1      , ONLY: mtopol,iat
      USE Zdaint_ag    , ONLY: daint1_ag,daint2_ag,D1_01,D1_bc,D1_09,aint_01b,aint_bc,aint_09b, &
                               daint1_ag_bc,aint_01d,aint_cm,aint_09d,daint2_ag_cm
      USE Zflowregime  , ONLY: alphag_bc,alphag_cm
      USE Ziat         , ONLY: ia,dbubble_init
      USE Ztimecon     , ONLY: alpha_min
      USE Ziat         , ONLY: s_bubble_diameter,r_ddrop,r_db_min,r_db_max,r_dh_hibiki
      USE Zvector      , ONLY: ug_o,ul_o
      USE Zturb        , ONLY: turb_dp      
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER i
      REAL(8) jg
      REAL(8) Dmin,lo,reb,dh
      REAL(8) deltarho,Dbub_max
      REAL(8) jg01,jgbc,jg09,gg     
!
!.....Topology selection with gamma_void, gradient of void fraction and alphag, void fraction
!
!
!.....1.Do not use topolgy map if mtopol=0. All regimes are 11.
!.....2.Use simple topology map if mtopol=1. All regimes defined as 11,12,13 by void fraction.
!
      IF(mtopol.eq.1) THEN
         DO i=1,ncell_fluid
            IF(cell%alphag(i).le.alphag_bc)THEN
               cell%regime(i)=11
            ELSEIF(cell%alphag(i).lt.alphag_cm)THEN
               cell%regime(i)=12
            ELSE
               cell%regime(i)=13
            ENDIF
         ENDDO
      ELSE
         cell%regime(:)=11
      ENDIF
!
!.....Calculate gamma_void, gradient of void fraction
!
      CALL initialize_topology(0)                         
!
!.....3. Use full topology map, if mtopol=2. All regimes are 11,21,3, 12,22,3, 13,23,3.
!      
      IF(mtopol.eq.2) CALL initialize_topology(1)     
!     
!.....Bubble Diameter Model Selection
!
! 
!.....'yoneda'          (default, Nuclear Engineering and Design 217 (2002) 267-281)
!
      IF(s_bubble_diameter.eq.'yoneda')THEN
         DO i=1,ncell_fluid
            cell%D1(i)= 10.06d0*(1.0e5/cell%p(i))**0.098d0  &
                       *SQRT(cell%sigma(i)/9.8d0/MAX(alpha_min,(cell%rhol(i)-cell%rhog(i))))
            D1_01(i)   =cell%D1(i)*0.1d0**0.35
            D1_bc(i)   =cell%D1(i)*0.118d0**0.35
            D1_09(i)   =D1_bc(i)
            cell%D1(i)=cell%D1(i)*(MIN(cell%alphag(i),0.118d0))**0.35d0
         ENDDO
!         
!.....'trac'            (TRAC-M code model)
!
      ELSEIF(s_bubble_diameter.eq.'trac')THEN
         DO i=1,ncell_fluid
            cell%D1(i)=2.0d0*SQRT(cell%sigma(i)/9.8d0/MAX(alpha_min,(cell%rhol(i)-cell%rhog(i))))
            D1_01(i)  =cell%D1(i)
            D1_bc(i)  =cell%D1(i)
            D1_09(i)  =cell%D1(i)
         ENDDO
!
!.....'modified_yoneda' (modified yoneda model for DOBO simulation by hkcho)
!
      ELSEIF(s_bubble_diameter.eq.'modified_yoneda')THEN
         DO i=1,ncell_fluid
            cell%D1(i)=10.06d0*(1.0e5/cell%p(i))**0.098d0*SQRT(cell%sigma(i)/9.8d0/ &
                       MAX1(alpha_min,(cell%rhol(i)-cell%rhog(i))))
            Dmin=2.0d0*SQRT(cell%sigma(i)/9.8d0/MAX1(alpha_min,(cell%rhol(i)-cell%rhog(i))))
            D1_01(i)  =Dmin
            D1_bc(i)  =cell%D1(i)*0.25d0**0.35d0
            D1_09(i)  =D1_bc(i)
            cell%D1(i)=cell%D1(i)*(DMIN1(cell%alphag(i),0.250))**0.35d0
            IF(cell%alphag(i).le.0.2d0)THEN
               cell%D1(i)=Dmin
            ELSEIF(cell%alphag(i).le.0.25d0)THEN
               cell%D1(i)=Dmin*(0.25d0-cell%alphag(i))/0.05d0+cell%D1(i)*(cell%alphag(i)-0.2d0)/0.05d0
            ELSE
               cell%D1(i)=cell%D1(i)
            ENDIF         
         ENDDO 
!
!.....'viswanathan'     (Canadian Journal of Chemical Engineering, 1985)
!
      ELSEIF(s_bubble_diameter.eq.'viswanathan')THEN            
         DO i=1,ncell_fluid
            jg=cell%alphag(i)*ug_o(i)
            jg01=0.1d0*ug_o(i)
            jgbc=alphag_bc*ug_o(i)
            jg09=0.9d0*ug_o(i)
            gg=cell%sigma(i)**3/cell%rhol(i)**3/9.81d0**2
            cell%D1(i)=LOG(MAX(0.000001d0,7.768d0*jg**0.26d0))  *(gg/MAX(0.000001d0,jg)**2)**0.2d0
            D1_01(i)  =LOG(MAX(0.000001d0,7.768d0*jg01**0.26d0))*(gg/MAX(0.000001d0,jg01)**2)**0.2d0
            D1_bc(i)  =LOG(MAX(0.000001d0,7.768d0*jgbc**0.26d0))*(gg/MAX(0.000001d0,jgbc)**2)**0.2d0
            D1_09(i)  =LOG(MAX(0.000001d0,7.768d0*jg09**0.26d0))*(gg/MAX(0.000001d0,jg09)**2)**0.2d0
         ENDDO 
!
!.....'hibiki'          (Chemical Engineering Science, 61, 7979-7990, 2006)
!
      ELSEIF(s_bubble_diameter.eq.'hibiki')THEN 
         dh=r_dh_hibiki
         CALL hibiki_bubble_diameter(cell%alphag,alphag_bc,    &
              dh,ul_o,ug_o,cell%rhol,cell%rhog,                &
              cell%lviscosl,cell%lviscosg,cell%sigma,          &
              cell%D1,D1_01,D1_bc,D1_09,ncell_fluid)
!
!.......'Yun'
!
      ELSEIF(s_bubble_diameter.eq.'yun')THEN  !for SUBO,DEBORA
         DO i=1,ncell_fluid
            lo=dsqrt(cell%sigma(i)/9.8d0/dmax1(alpha_min,(cell%rhol(i)-cell%rhog(i))))
            Reb=turb_dp(i)**(1.0d0/3.0d0)*lo**(4.0d0/3.0d0)*cell%rhol(i)/cell%lviscosl(i)
!            Reb=1.d-3**(1.0d0/3.0d0)*lo**(4.0d0/3.0d0)*cell(i)%rhol/cell(i)%lviscosl
            IF(Reb.eq.0.0d0)Then
               cell%D1(i)=0.0d0
            ELSE
               cell%D1(i)=32.39d0*lo*Reb**(-0.529d0)*(cell%rhol(i)/cell%rhog(i))**(0.06d0)
               D1_01(i)  =cell%D1(i)*0.1d0**0.289d0
               D1_bc(i)  =cell%D1(i)*alphag_bc**0.289d0
               D1_09(i)  =cell%D1(i)*0.9d0**0.289d0
               cell%D1(i)=cell%D1(i)*cell%alphag(i)**0.289d0
            ENDIF   
          ENDDO    
      ELSEIF(s_bubble_diameter.eq.'yun_old')THEN  !for core_catcher
         DO i=1,ncell_fluid
            lo=sqrt(cell%sigma(i)/9.8d0/max(alpha_min,(cell%rhol(i)-cell%rhog(i))))
            Reb=turb_dp(i)**(1.0d0/3.0d0)*lo**(4.0d0/3.0d0)*cell%rhol(i)/cell%lviscosl(i)
!            Reb=1.d-3**(1.0d0/3.0d0)*lo**(4.0d0/3.0d0)*cell(i)%rhol/cell(i)%lviscosl
            IF(Reb.eq.0.0d0)Then
               cell%D1(i)=0.0d0
            ELSE
               cell%D1(i)=39.32d0*lo*Reb**(-0.696d0)*(cell%rhol(i)/cell%rhog(i))**(0.571d0)
               D1_01(i)  =cell%D1(i)*0.1d0**0.36d0
               D1_bc(i)  =cell%D1(i)*alphag_bc**0.36d0
               D1_09(i)  =cell%D1(i)*0.9d0**0.36d0
               cell%D1(i)=cell%D1(i)*cell%alphag(i)**0.36d0
            ENDIF   
          ENDDO                 
!                      
      ENDIF           
      
!     
!.....Calculate interfacial area (IA) 
!
!
!........Bubble diameters and gradients
!
      DO i=1,ncell_fluid
         cell%D1(i)=MIN(r_db_max,MAX(r_db_min,cell%D1(i)))
         D1_01(i)  =MIN(r_db_max,MAX(r_db_min,D1_01(i)))
         D1_bc(i)  =MIN(r_db_max,MAX(r_db_min,D1_bc(i)))
         D1_09(i)  =MIN(r_db_max,MAX(r_db_min,D1_09(i)))       
         daint1_ag(i)=6.0d0/cell%D1(i)
         cell%aint1(i)=daint1_ag(i)*MAX(1.0d-8,cell%alphag(i))
         cell%aint1(i)=MIN(cell%aint1(i),6.0d0*alphag_bc/D1_bc(i))
      ENDDO               
!         
      DO i=1,ncell_fluid
         aint_01b(i)    =6.0d0/D1_01(i)*0.1d0
         aint_bc(i)     =6.0d0/D1_bc(i)*alphag_bc
         aint_09b(i)    =6.0d0/D1_09(i)*alphag_bc
         daint1_ag_bc(i)=6.0d0/D1_bc(i)
         IF(cell%alphag(i).gt.alphag_bc) daint1_ag(i)=0.0d0
         IF(IAT.eq.0)THEN
            cell%D1(i)=dbubble_init
            cell%aint1(i)=6.0d0*cell%alphag(i)/cell%D1(i)
         ENDIF            
         IF(IAT.gt.0)THEN
            cell%aint1(i)=ia(i)
!
!........Set minimum & maximum bubble size
!
            deltarho  =MAX(1.0d-5,cell%rhol(i)-cell%rhog(i)) 
            Dbub_max  =40.d0*SQRT(cell%sigma(i)/9.806d0/deltarho)
            cell%D1(i)=MAX(r_db_min,6.0d0*cell%alphag(i)/cell%aint1(i))
            cell%D1(i)=MIN(Dbub_max,cell%D1(i))  
!            
            cell%aint1(i)=6.0d0*cell%alphag(i)/cell%D1(i)
         ENDIF
         ia(i)=cell%aint1(i)
!
!........Droplet diameters and gradients
!
         cell%D2(i)=r_ddrop
         daint2_ag(i)=6.0d0/cell%D2(i)
         cell%aint2(i)=daint2_ag(i)*MAX(1.0d-8,1.0d0-cell%alphag(i))
         cell%aint2(i)=MIN(cell%aint2(i),6.0d0*(1.0d0-alphag_cm)/cell%D2(i))
         daint2_ag(i)=-daint2_ag(i)
         aint_01d(i)=6.0d0/cell%D2(i)
         daint2_ag_cm(i)=-aint_01d(i)
         aint_cm(i)=aint_01d(i)*(1.0d0-alphag_cm)
         aint_09d(i)=aint_01d(i)*0.1d0
         aint_01d(i)=aint_01d(i)*0.9d0
      ENDDO
!
      cell%entr(:)=0.0d0
      cell%dentr(:)=0.0d0
      cell%yeta(:)=0.0d0      
!
      END SUBROUTINE int_area
