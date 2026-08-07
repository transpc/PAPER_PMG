!
      SUBROUTINE ncg_cell 
!
!     This routine noncondensible gas properties for each computing cells      
!
      USE Zmpi         , ONLY: ncell_fp
      USE Zzone        , ONLY: ncell_fluid
      USE Zparam       , ONLY: nb_max
      USE Zrv_model    , ONLY: rv_ht_w
      USE Zconst1      , ONLY: wconden,rv_htmodel_forCFD
      USE STM_TBL_cupid  , ONLY: wmole,dcvax,cvaox,uaox,rmolg,advn
      USE Zbc_index    , ONLY: nvin,npin
      USE Zncg         , ONLY: n_ncg_sp,ncg_species,tao,qn_cell_o,qn_cell0,           &
                               cvao_cell,uao_cell,dcva_cell,ra_cell,qn_cell,advn_cell,&
                               cvao_nvin,uao_nvin,dcva_nvin,ra_nvin,qn_nvin,          &
                               cvao_npin,uao_npin,dcva_npin,ra_npin,qn_npin,wmole_gas
      USE Zncg, ONLY: i_ncg_vis
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i,n,inc,ncg
      LOGICAL,SAVE :: initial=.true.
      REAL(8) :: t1,t2,t3,t4,t5
      REAL(8) :: xmsum
      REAL(8) :: cvaox0(n_ncg_sp),uaox0(n_ncg_sp),dcvax0(n_ncg_sp),wmole0(n_ncg_sp),advn0(n_ncg_sp)
!
!     wmolea  molecular mass of non-condensible gas.
!     rax     gas constant of non-condensible gas.
!     dcvax   same as above.
!     cvaox   same as above.
!     uaox    term in u = uao + integral (cv*dt) where u is internal
!             energy.
!     tao     term in cv = cvao + dcva*(t - tao) where cv is heat capacity
!             and t is temperature.
!     dconst  DIFfusion coefficient at reference conditions for non-
!             condensible gasses and steam.
!     noncn   number of non-condensible gasses.
!     prop    array for sth2x CALLs, also USEd for scratch.
!     s       same as above.
!
!
!     Data Statement:  Constant for evaluation of the DIFfusion            
!                      Coefficient of NC gas in Water Vapor.               
!                                                                         
!                      diffc = Dconst * T**1.75 / P                        
!                                                                       
!     Ref:  eq. 11-4.1 of "Properties of Gases and Liquids"           
!           by Reid, Praudnitz & Sherwood.                            
!           3rd edition, McGraw-Hill Book Co., 1977.                  
!                                                                       
!     NC gases are: Helium, Hydrogen, Nitrogen, Krypton, Xenon, Air, Argon, SF6.                                    
!
      IF(initial)THEN
         n=ncell_fluid
         ALLOCATE(cvao_cell(n),uao_cell(n),dcva_cell(n),ra_cell(n),advn_cell(n))
         ALLOCATE(qn_cell(ncell_fp,n_ncg_sp),qn_cell_o(ncell_fp,n_ncg_sp))
         IF(i_ncg_vis.ge.1.or.ABS(wconden).eq.1.or.ABS(wconden).eq.2.or.rv_ht_w.eq.1.or.rv_htmodel_forCFD.gt.0)ALLOCATE(wmole_gas(n))
         n=nb_max
         ALLOCATE(cvao_nvin(n),uao_nvin(n),dcva_nvin(n),ra_nvin(n))
         ALLOCATE(cvao_npin(n),uao_npin(n),dcva_npin(n),ra_npin(n))
         initial=.false.
!
         DO i=1,n_ncg_sp
            qn_cell(:,i)=qn_cell0(i)
         ENDDO
!
      ENDIF
!
      tao=250.0d0 
!
!DIR$ NOVECTOR
      DO inc=1,n_ncg_sp
         ncg=ncg_species(inc)
         cvaox0(inc)=cvaox(ncg)
         uaox0(inc)=uaox(ncg)
         dcvax0(inc)=dcvax(ncg)
         wmole0(inc)=wmole(ncg)
         advn0(inc)=advn(ncg)
      END DO
!
!DIR$ SIMD
      DO i=1,ncell_fluid
         t1=0.d0
         t2=0.d0
         t3=0.d0
         t4=0.d0
         t5=0.d0
         DO inc=1,n_ncg_sp
            t1=t1+cvaox0(inc)*qn_cell(i,inc)
            t2=t2+uaox0(inc) *qn_cell(i,inc)
            t3=t3+dcvax0(inc)*qn_cell(i,inc)
            t4=t4+rmolg/wmole0(inc)*qn_cell(i,inc)
            t5=t5+advn0(inc)*qn_cell(i,inc)
         ENDDO
         cvao_cell(i)=t1
         uao_cell(i) =t2
         dcva_cell(i)=t3
         ra_cell(i)  =t4
         advn_cell(i)=t5
      END DO
!
      DO i=1,nvin
         t1=0.d0
         t2=0.d0
         t3=0.d0
         t4=0.d0
!DIR$ NOVECTOR
         DO inc=1,n_ncg_sp
            t1=t1+cvaox0(inc)*qn_nvin(i,inc)
            t2=t2+uaox0(inc) *qn_nvin(i,inc)
            t3=t3+dcvax0(inc)*qn_nvin(i,inc)
            t4=t4+rmolg/wmole0(inc)*qn_nvin(i,inc)
         ENDDO
         cvao_nvin(i)=t1
         uao_nvin(i) =t2
         dcva_nvin(i)=t3
         ra_nvin(i)  =t4
      END DO
! 
      DO i=1,npin
         t1=0.d0
         t2=0.d0
         t3=0.d0
         t4=0.d0
!DIR$ NOVECTOR
         DO inc=1,n_ncg_sp
            t1=t1+cvaox0(inc)*qn_npin(i,inc)
            t2=t2+uaox0(inc) *qn_npin(i,inc)
            t3=t3+dcvax0(inc)*qn_npin(i,inc)
            t4=t4+rmolg/wmole0(inc)*qn_npin(i,inc)
         ENDDO
         cvao_npin(i)=t1
         uao_npin(i) =t2
         dcva_npin(i)=t3
         ra_npin(i)  =t4
      END DO
!
!.....Ave. molar weight of NCG
!      
      IF(i_ncg_vis.ge.1.or.ABS(wconden.eq.1).or.ABS(wconden.eq.2).or.rv_ht_w.gt.0.or.rv_htmodel_forCFD.gt.0)THEN
         DO i=1,ncell_fluid 
            xmsum=0.d0
            DO inc=1,n_ncg_sp    
               xmsum=xmsum+qn_cell(i,inc)/wmole0(inc)
            ENDDO
            wmole_gas(i)=1.d0/xmsum
         ENDDO
      ENDIF
!
      END SUBROUTINE ncg_cell 




