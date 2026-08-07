!
      SUBROUTINE udfn_porous_property
!
!     Define Aporous (only when "udfl_porous_property" is used)
!    
      USE Zzone        , ONLY: ncell_fluid
      USE Zcore        , ONLY: myrank,np
      USE Zparam       , ONLY: pi
      USE Zconst1      , ONLY: vv_prob
      USE Zcoord3      , ONLY: aporous,porosity,vol
      USE Zvec_geo     , ONLY: sa_nf
      USE Zio_unit     , ONLY: unit_log
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER i
      REAL(8) corevols,aporoussums      
      REAL(8) porous_s,nrod_tot,aporous_tot !mcc-pik
      LOGICAL,SAVE::initial   
      DATA initial/.true./  
!
      nrod_tot=0.0d0
      aporous_tot=0.0d0  
!     
!.....Aporous [m2] is the actural area of a solid cell in porous cases
!         
      IF(vv_prob.eq.'PAFS-POOL')THEN
         DO i=1,ncell_fluid
            Aporous(i)=0.032d0
         ENDDO    
      ENDIF  
!
      IF(vv_prob.eq.'stern') THEN
         DO i=1,ncell_fluid
            Aporous(i)=20.58d0 
         ENDDO    
      ENDIF  
!
      IF(vv_prob.eq.'annul_porous') then
      ENDIF
!      
      IF(vv_prob.eq.'block_porous') then
      ENDIF
!         
      IF(vv_prob.eq.'apr1400_lbloca') then
         CALL max_nf(sa_nf,aporous)
      ENDIF
!
      IF(vv_prob.eq.'copain_porous')THEN
         DO i=1,ncell_fluid
            aporous(i)=0.005d0
         ENDDO    
      ENDIF
!      
!.....ATLAS, use geometrical data obtained from Dr. B.U.Bae
!
      IF(vv_prob.eq.'atlas_mc_porous'.or.vv_prob.eq.'pwr_mc_poro'.or.vv_prob.eq.'apr1400_mc_poro'.or.vv_prob.eq.'opr1000_mc_poro')THEN
          IF(initial)then
             initial=.FALSE.
             corevols=0.0d0 
                 DO i=1,ncell_fluid
                    porous_s=1.0d0-porosity(i)
                    corevols=corevols+vol(i)*porous_s
                 ENDDO  
             IF(np.gt.1) CALL allreducei_r1(corevols)
             aporous(:)=0.0d0
             aporoussums=0.0d0
             IF(vv_prob.eq.'atlas_mc_porous')then
                 DO i=1,ncell_fluid
                    porous_s=1.0d0-porosity(i)
                    aporous(i)=2.0d0*pi*(0.095d0/2.d0)*1.9d0*390.0d0*vol(i)*porous_s/corevols !must be saved
                    aporoussums=aporoussums+aporous(i)
                 ENDDO 
              ELSEIF(vv_prob.eq.'pwr_mc_poro'.or.vv_prob.eq.'apr1400_mc_poro'.or.vv_prob.eq.'opr1000_mc_poro')THEN
                 DO i=1,ncell_fluid
                    porous_s=1.0d0-porosity(i)
                    aporous(i)=2.0d0*pi*(0.095d0/2.d0)*3.8d0*241*236*vol(i)*porous_s/corevols !must be saved
                    aporoussums=aporoussums+aporous(i)
                 ENDDO                    
              ENDIF    
             IF(np.gt.1) CALL allreducei_r1(aporoussums)
             IF(myrank.eq.0)then
                WRITE(*,"(11x,a,1pe17.5,a)")'core_solid_vol=',corevols,'m3'   
                WRITE(*,"(11x,a,1pe17.5,a)")'aporous_sum=',aporoussums,'m2'             
                WRITE(unit_log,"(11x,a,1pe17.5,a)")'core_solid_vol=',corevols,'m3'   
                WRITE(unit_log,"(11x,a,1pe17.5,a)")'aporous_sum=',aporoussums,'m2' 
             ENDIF   
          ELSE
          ENDIF 
      ENDIF
!       
      END SUBROUTINE udfn_porous_property
