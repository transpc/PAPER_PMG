!
      SUBROUTINE udfn_subchannel_CHF
!
      USE Zinterface
      USE VOL_DATA       , ONLY: cell
      USE Zmpi           , ONLY: ncell_fp
      USE Zzone          , ONLY: ncell_fluid
      USE Zcore          , ONLY: np,myrank
      USE Znum_cell      , ONLY: i_neigh,neigh
      USE Ztimecon       , ONLY: time 
      USE Zparam         , ONLY: ns,pi
      USE Zporous    , ONLY: chn_type
      USE Zconst1        , ONLY: vv_prob
      USE Zconst2        , ONLY: hydraulicd
      USE Zvector        , ONLY: vl_n,vg_n
      USE Zvec_geo       , ONLY: xn_nf
      ! CUPID-RV
      USE zrv_ncell      , ONLY: ncell_fuel_rod,p3d_cupid,cupid_cell_hts2d,dnbr_ce1,dnbr_ce2,master_to_rod,dnbr_cupid1,dnbr_cupid2
      USE Zrv_hts_2d     , ONLY: ri_2d,nr_2d
      USE Zcoord2        , ONLY: cell_leng
      USE Zrv_model      , ONLY: rv_model
      USE MASTER4        , ONLY: NZ_TH,NXYF,NPINX,PIN3D_TH
      USE Zrv_mpi        , ONLY: jperm_fuel_rod
      USE Zwall_HTC      , ONLY: qflux_l0    
!
      IMPLICIT NONE
!.....Local variables
      INTEGER i,j,j0
      REAL(8) b1, b2, b3, b4, b5,b6, b7,b8 
      REAL(8) hd_matrix 
      ! CUPID-RV
      INTEGER m1,m2,m4
      INTEGER ij,im,ii,iz,ic,irod
      LOGICAL,SAVE :: init_dnbr=.true.
      REAL(8) dm1
      REAL(8) rhog,alphag,hgsat,hg,vg
      REAL(8) rhol,alphal,hlsat,hl,vl
      REAL(8) quals,p
      REAL(8) min_dnbr1
      REAL(8) min_dnbr2
!.....Local arrays
      REAL(8) xn(ns),yn(ns)
      REAL(8) dm(ncell_fp),dm0(ncell_fp)
      REAL(8) rhog0(ncell_fp),alphag0(ncell_fp),hgsat0(ncell_fp),hg0(ncell_fp),vg0(ncell_fp)
      REAL(8) rhol0(ncell_fp),alphal0(ncell_fp),hlsat0(ncell_fp),hl0(ncell_fp),vl0(ncell_fp)
      REAL(8) p0(ncell_fp),quals0(ncell_fp)
!
      REAL(8) Aprime(ncell_fuel_rod)
      REAL(8) A(ncell_fuel_rod), B(ncell_fuel_rod),C(ncell_fuel_rod),G(ncell_fuel_rod)
      REAL(8) q_CHF_CE(ncell_fuel_rod), Hfg(ncell_fuel_rod), Qual(ncell_fuel_rod),hmix(ncell_fuel_rod),q_cupid(ncell_fuel_rod)

!
      b1= 2.8922d-3
      b2=-0.50749d0
      b3= 405.320d0
      b4=-9.9290d-2
      b5=-0.67757d0
      b6= 6.8235d-4
      b7= 3.1240d-4
      b8=-8.3245d-2
!      dm=heated area만 고려한 hydraulic diameter

      !IF(vv_prob.eq.'single_assem' .or. vv_prob.eq.'APR1400_fullcore' .or. vv_prob.eq.'OPR1000_fullcore')THEN
      IF(vv_prob.eq.'OPR1000_fullcore_modmesh02_rv' .or. &
         vv_prob.eq.'OPR1000_fullcore_modmesh02_rv_vessel')THEN
!        CE-1 Correlation

         hd_matrix=0.012637d0  !heated diameter of center subchannel=equivalent heated diameter of matrix subchannel
         IF(vv_prob.eq.'OPR1000_fullcore'             .or. &
            vv_prob.eq.'OPR1000_quarter_core'         .or. &
            vv_prob.eq.'OPR1000_fullcore_modmesh02_rv'.or. &
            vv_prob.eq.'OPR1000_fullcore_modmesh02_rv_vessel') hd_matrix=0.01264d0

         DO i=1,ncell_fluid
            dm(i)=0.0d0
            IF(vv_prob.eq.'APR1400_fullcore')THEN   !180705 KSB modi
               IF(chn_type(i).eq.2)THEN !side
                   dm(i)=0.012637d0
               ELSEIF(chn_type(i).eq.3)THEN !corner
                   dm(i)=0.012637d0
               ELSEIF(chn_type(i).eq.5)THEN !guide side
                   dm(i)=0.014558d0
               ELSEIF(chn_type(i).eq.6)THEN !guide corner
                   dm(i)=0.016426d0
               ELSEIF(chn_type(i).gt.6)THEN !water gap, shroud
                   dm(i)=1.d0
               ELSE
                   dm(i)=hydraulicd(i)
               ENDIF
            ELSEIF(vv_prob.eq.'OPR1000_fullcore')THEN
               IF(chn_type(i).eq.2)THEN !side
                   dm(i)=0.01264d0
               ELSEIF(chn_type(i).eq.3)THEN !corner
                   dm(i)=0.01264d0
               ELSEIF(chn_type(i).eq.5)THEN !guide side
                   dm(i)=0.013855d0
               ELSEIF(chn_type(i).eq.6)THEN !guide corner
                   dm(i)=0.016922d0
               ELSEIF(chn_type(i).gt.6)THEN !water gap, shroud
                   dm(i)=1.d0
               ELSE
                   dm(i)=hydraulicd(i)
               ENDIF
            ELSEIF(vv_prob.eq.'OPR1000_fullcore_modmesh02_rv' .or. &
                   vv_prob.eq.'OPR1000_fullcore_modmesh02_rv_vessel')THEN
               IF(chn_type(i).eq.2)THEN !side
                   dm(i)=0.01264d0
               ELSEIF(chn_type(i).eq.3)THEN !corner
                   dm(i)=0.01264d0
               ELSEIF(chn_type(i).eq.5)THEN !guide side
                   dm(i)=0.013855d0
               ELSEIF(chn_type(i).eq.6)THEN !guide corner
                   dm(i)=0.016922d0
               ELSEIF(chn_type(i).gt.6)THEN !water gap, shroud
                   dm(i)=1.d0
               ELSE
                   dm(i)=hydraulicd(i)
               ENDIF
            ENDIF
         ENDDO

         IF(np.gt.1) THEN
            CALL communicate_1d(dm         , &
                                cell%rhol  , &
                                cell%rhog  , &
                                cell%alphal, &
                                cell%alphag, &
                                cell%p     , &
                                cell%hgsat , &
                                cell%hlsat  ) 
            CALL communicate_1d(cell%quals , &
                                cell%hl    , &
                                   cell%hg    ) 
            CALL communicate_2d(vl_n, &
                                vg_n)    
         ENDIF
                
         DO m1=1,ncell_fluid
            IF(chn_type(m1).ne.0)then
               m2=0
!..............Get all the xn(j,1) for cell m1
               CALL  get_scalar_variable_n_i_ndim(xn_nf,xn,m1,1)
               j0=i_neigh(m1)-1
               DO j=i_neigh(m1),i_neigh(m1+1)-1
                  IF(xn(j-j0).gt.0.5d0) m2=neigh(j)
               ENDDO
               IF(m2.ne.0)then
                  dm0(m1)    =(dm(m1)          +dm(m2)         )*0.5d0
                  rhog0(m1)  =(cell%rhog(m1)   +cell%rhog(m2)  )*0.5d0
                  rhol0(m1)  =(cell%rhol(m1)   +cell%rhol(m2)  )*0.5d0
                  alphag0(m1)=(cell%alphag(m1) +cell%alphag(m2))*0.5d0
                  alphal0(m1)=(cell%alphal(m1) +cell%alphal(m2))*0.5d0
                  p0(m1)     =(cell%p(m1)      +cell%p(m2)     )*0.5d0
                  hgsat0(m1) =(cell%hgsat(m1)  +cell%hgsat(m2) )*0.5d0
                  hlsat0(m1) =(cell%hlsat(m1)  +cell%hlsat(m2) )*0.5d0
                  quals0(m1) =(cell%quals(m1)  +cell%quals(m2) )*0.5d0
                  hg0(m1)    =(cell%hg(m1)     +cell%hg(m2)    )*0.5d0
                  hl0(m1)    =(cell%hl(m1)     +cell%hl(m2)    )*0.5d0
                  vl0(m1)    =(vl_n(m1,3)      +vl_n(m2,3)     )*0.5d0
                  vg0(m1)    =(vg_n(m1,3)      +vg_n(m2,3)     )*0.5d0
               ELSE             
                  dm0(m1)    =dm(m1)                       
                  rhog0(m1)  =cell%rhog(m1)  
                  rhol0(m1)  =cell%rhol(m1)  
                  alphag0(m1)=cell%alphag(m1)
                  alphal0(m1)=cell%alphal(m1)
                  p0(m1)     =cell%p(m1)     
                  hgsat0(m1) =cell%hgsat(m1) 
                  hlsat0(m1) =cell%hlsat(m1) 
                  quals0(m1) =cell%quals(m1) 
                  hg0(m1)    =cell%hg(m1)    
                  hl0(m1)    =cell%hl(m1)    
                  vl0(m1)    =vl_n(m1,3)     
                  vg0(m1)    =vg_n(m1,3)     
               ENDIF   
            ELSE
               dm0(m1)    =dm(m1)
               rhog0(m1)  =cell%rhog(m1)
               rhol0(m1)  =cell%rhol(m1)
               alphag0(m1)=cell%alphag(m1)
               alphal0(m1)=cell%alphal(m1)
               p0(m1)     =cell%p(m1)
               hgsat0(m1) =cell%hgsat(m1)
               hlsat0(m1) =cell%hlsat(m1)
               quals0(m1) =cell%quals(m1)
               hg0(m1)    =cell%hg(m1)
               hl0(m1)    =cell%hl(m1)
               vl0(m1)    =vl_n(m1,3)
               vg0(m1)    =vg_n(m1,3)
            ENDIF
         ENDDO

         IF(np.gt.1) THEN
            CALL communicate_1d(dm0    , &
                                rhog0  , &
                                rhol0  , &
                                alphal0, &
                                alphag0, &
                                p0     , &
                                hgsat0 , &
                                hlsat0    ) 
            CALL communicate_1d(quals0 , &
                                hl0    , &
                                hg0    , &
                                vl0    , &
                                vg0       ) 
         ENDIF

         IF(init_dnbr)then
            ALLOCATE(dnbr_ce1(ncell_fuel_rod))
            ALLOCATE(dnbr_ce2(ncell_fuel_rod))
            dnbr_ce1=0.0d0
            dnbr_ce2=0.0d0
            init_dnbr=.false.
         ENDIF

if(0)then
         IF(rv_model.eq.0)then
         DO iz=1,NZ_TH
            DO im=1,NXYF
               DO ij=1,NPINX
                  DO ii=1,NPINX
                     ic=master_to_rod(ii,ij,im,iz)
                     DO irod=1,ncell_fuel_rod
                        i=jperm_fuel_rod(irod)
                        IF(i.eq.ic)then
                           p3d_cupid(irod)=PIN3D_TH(ii,ij,im,iz)
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         ENDIF
endif

         DO i=1,ncell_fuel_rod
            dnbr_ce1(i)=0.0d0
            m1=cupid_cell_hts2d(i)
            m4=0
!...........Get all the xn(j,1),xn(j,2) for cell m1
            CALL get_scalar_variable_n_i_ndim(xn_nf,yn,m1,2)
            j0=i_neigh(m1)-1
            DO j=i_neigh(m1),i_neigh(m1+1)-1
               IF(yn(j-j0).gt.0.5d0) m4=neigh(j)
            ENDDO
            IF(m4.ne.0)then
               dm1   =(dm0(m1)     +dm0(m4)    )*0.5d0
               rhog  =(rhog0(m1)   +rhog0(m4)  )*0.5d0
               rhol  =(rhol0(m1)   +rhol0(m4)  )*0.5d0
               alphag=(alphag0(m1) +alphag0(m4))*0.5d0
               alphal=(alphal0(m1) +alphal0(m4))*0.5d0
               p     =(p0(m1)      +p0(m4)     )*0.5d0
               hgsat =(hgsat0(m1)  +hgsat0(m4) )*0.5d0
               hlsat =(hlsat0(m1)  +hlsat0(m4) )*0.5d0
               quals =(quals0(m1)  +quals0(m4) )*0.5d0
               hg    =(hg0(m1)     +hg0(m4)    )*0.5d0
               hl    =(hl0(m1)     +hl0(m4)    )*0.5d0
               vl    =(vl0(m1)     +vl0(m4)    )*0.5d0
               vg    =(vg0(m1)     +vg0(m4)    )*0.5d0     
            ELSE
               dm1   =dm0(m1)     
               rhog  =rhog0(m1)   
               rhol  =rhol0(m1)   
               alphag=alphag0(m1) 
               alphal=alphal0(m1) 
               p     =p0(m1)      
               hgsat =hgsat0(m1)  
               hlsat =hlsat0(m1)  
               quals =quals0(m1)  
               hg    =hg0(m1)     
               hl    =hl0(m1)     
               vl    =vl0(m1)     
               vg    =vg0(m1)     
            ENDIF

if(0)then
               dm1   =dm0(m1)
               rhog  =cell%rhog(m1)
               rhol  =cell%rhol(m1)
               alphag=cell%alphag(m1)
               alphal=cell%alphal(m1)
               p     =cell%p(m1)
               hgsat =cell%hgsat(m1)
               hlsat =cell%hlsat(m1)
               quals =cell%quals(m1)
               hg    =cell%hg(m1)
               hl    =cell%hl(m1)
               vl    =vl_n(m1,3)
               vg    =vg_n(m1,3)
endif

            G(i)=(rhol*vl*alphal+rhog*vg*alphag)/1356.229d0 !British unit(lb/hr*ft^2)
            !Aprime(i)=b1*(hydraulicd(i)/dm(i))**b2     ! original
            Aprime(i)=b1*(dm1/hd_matrix)**b2      ! 180709 KSB modi
            A(i)=(b3+b4*p/6894.8d0)*G(i)**(b5+b6*p/6894.8d0)
            B(i)=G(i)*(hgsat-hlsat)/2326.011d0 !latent heat[Btu/lb]
            C(i)=G(i)**(b7*p/6894.8d0+b8*G(i))
            hmix(i)=hl*(1.d0-quals)+hg*quals
            Qual(i)=(hmix(i)-hlsat)/(hgsat-hlsat)
            q_CHF_CE(i)=Aprime(i)*(A(i)-B(i)*Qual(i))/C(i)*3.1546d6
            Hfg(i)=hgsat-hlsat
            
            IF(p3d_cupid(i).gt.1.e-3)then
               IF(rv_model.eq.0) then
                  q_cupid(i)=p3d_cupid(i)/(2.0*pi*ri_2d(nr_2d)*cell_leng(m1,3))
                  dnbr_ce1(i)=q_CHF_CE(i)/q_cupid(i)
                  dnbr_cupid1(m1)=dnbr_ce1(i)
               ENDIF
               IF(rv_model.eq.1) then
                  !cupid-rv
                  dnbr_ce1(i)=q_CHF_CE(i)/(qflux_l0(i))
                  dnbr_cupid1(m1)=dnbr_ce1(i)

                  !master
                  q_cupid(i)=p3d_cupid(i)/(2.0*pi*ri_2d(nr_2d)*cell_leng(m1,3))
                  dnbr_ce2(i)=q_CHF_CE(i)/q_cupid(i)
                  dnbr_cupid2(m1)=dnbr_ce2(i)
               ENDIF
            ELSE
               dnbr_ce1(i)=10.d0
               dnbr_ce2(i)=10.d0
               dnbr_cupid1(m1)=10.d0
               dnbr_cupid2(m1)=10.d0
            ENDIF
         ENDDO

         min_dnbr1=1.e20
         min_dnbr2=1.e20
         DO i=1,ncell_fuel_rod
            IF(dnbr_ce1(i).lt.min_dnbr1)then
               min_dnbr1=dnbr_ce1(i)
            ENDIF
            IF(dnbr_ce2(i).lt.min_dnbr2)then
               min_dnbr2=dnbr_ce2(i)
            ENDIF
         ENDDO

         IF(np.gt.1) THEN
            CALL allreducei_min_r1(min_dnbr1)
            CALL allreducei_min_r1(min_dnbr2)
         ENDIF
         
         IF(myrank.eq.0)then
            write(5801,*)time,min_dnbr1,min_dnbr2
         ENDIF

      ENDIF
!     
      END SUBROUTINE udfn_subchannel_CHF
!
