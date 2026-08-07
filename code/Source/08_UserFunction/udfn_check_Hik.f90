!
      SUBROUTINE udfn_check_Hik
!
!     Plot Hil and Hig as a function of void fraction and temperature
!
      USE VOL_DATA                    
      USE Zmpi            , ONLY: ncell_fp
      USE Zzone           , ONLY: ncell_fluid
      USE Zcore           , ONLY: np
      USE Zparam          , ONLY: ndim
      USE Zdhda           , ONLY: dHldag,dHgdag,dHldtl,dHgdtg
      USE Zqvol           , ONLY: h_il,h_ig
      USE Ziat            , ONLY: r_dh_hibiki
      USE Zvector         , ONLY: vrel_o
      USE Zvector         , ONLY: vl_o,vg_o
      USE Zvoid           , ONLY: gamma_void
      USE Zrv_model       , ONLY: free_model,rv_ht_i
      USE Zio_unit        , ONLY: unit_log
!
      IMPLICIT NONE
!
      INTEGER i,j,ii,nag,nt,nt1
!      
      REAL(8) elsat,egsat,del,deg,dag,ddag
      REAL(8), ALLOCATABLE::Cdl(:),Cdd(:)
      REAL(8), ALLOCATABLE::drag_coeff_c1(:)
!
!.....Output files to make a x-y graph of Hik w.r.t. void fraction for different number (nt) of temperatures
!
      OPEN(800,FILE='Hil_ag.dat')
      OPEN(900,FILE='Hig_ag.dat')
!
!.....Output files to make a x-y graph of Hik w.r.t. temperature when alphag=0.0 and alphag=1.0
!
      OPEN(801,FILE='Hil_T_ag0_ag1.dat')
      OPEN(901,FILE='Hig_T_ag0_ag1.dat')
!
!.....Outout files to plot Hik surface w.r.t void fraction and temperature
!
      OPEN(802,FILE='Hil_T_ag.dat')
      OPEN(902,FILE='Hig_T_ag.dat')
!
!.....Outout files to plot Hik surface w.r.t void fraction and void fraction gradient
!
      OPEN(803,FILE='Hil_ag_dag.dat')
      OPEN(903,FILE='Hig_ag_dag.dat')
!
!.....Outout files to plot derivative of Hik surface w.r.t void fraction and temperature
!
      OPEN(804,FILE='dHlda_ag_t.dat')
      OPEN(904,FILE='dHgda_ag_t.dat')
!
!.....Outout files to plot derivative of Hik surface w.r.t void fraction and temperature
!
      OPEN(805,FILE='dHldt_ag_t.dat')
      OPEN(905,FILE='dHgdt_ag_t.dat')
!
!.....Outout files to plot Cd surface w.r.t void fraction and void fraction gradient
!
      OPEN(810,FILE='Cd_ag_dag.dat')
!
      IF(np.gt.1) STOP '## check_Hik works only for single CPU. ##'
!
      elsat=cell%elsat(1)
      egsat=cell%egsat(1)
!
!.....nag: Number of plots for void fraction
!.....nt:  Number of plots for temperature
!
      nag=51
      nt=21
      nt1=nt/2+1
!
      ALLOCATE(Cdl(nag*nt),Cdd(nag*nt))
      ALLOCATE(drag_coeff_c1(nag*nt))
!
!.....Increment of void fraction and temperature(internal energy)
!
      dag=1.0d0/(nag-1)
      del=80000.d0/(nt-1)
      deg=90000.d0/(nt-1)
      ddag=1.0d0/(nt-1)
!
      IF(ncell_fluid.lt.nt*nag)THEN
         WRITE(*,*) '## More than',nt*nag,' cells are required for check_Hik. ##'
         WRITE(unit_log,*) '## More than',nt*nag,' cells are required for check_Hik. ##'
         STOP
      ENDIF
!
!.....Set plotting points of temperature
!
      ii=0
      DO i=1,nt
         ii=1+(i-1)*nag
         cell%el(ii:ii+nag-1)=elsat-(nt1-i)*del
         cell%eg(ii:ii+nag-1)=egsat-(nt1-i)*deg
         gamma_void(ii:ii+nag-1)=ddag*(i-1)
      ENDDO
!
!.....Set plotting points of void fraction
!
      ii=0
      DO j=1,nt
         DO i=1,nag
            ii=ii+1
            cell%alphag(ii)=dag*(i-1)
            IF(cell%alphag(ii).eq.0.0d0) cell%alphag(ii)=1.5d-8
            IF(cell%alphag(ii).eq.1.0d0) cell%alphag(ii)=1.0d0-1.5d-8
            cell%alphal(ii)=1.0d0-cell%alphag(ii)
         ENDDO
      ENDDO
!
      IF(ndim.eq.2) THEN
         DO i=1,ncell_fp
            vl_o(i,1)=1.0d0
            vl_o(i,2)=1.0d0
            vg_o(i,1)=1.0d0
            vg_o(i,2)=1.0d0
         ENDDO
      ELSE
         DO i=1,ncell_fp
            vl_o(i,1)=1.0d0
            vl_o(i,2)=1.0d0
            vl_o(i,3)=1.0d0
            vg_o(i,1)=1.0d0
            vg_o(i,2)=1.0d0
            vg_o(i,3)=1.0d0
         ENDDO
      ENDIF
      vrel_o(:)=1.0d0
      r_dh_hibiki=0.1d0
!
!.....Call the related subroutines for the caculation of Hik
!
      CALL property_calc(1)
      CALL int_area
      IF(free_model)CALL int_htc
      CALL int_swap(2)
      IF(rv_ht_i.gt.0)CALL rv_int_ht
      CALL int_swap(22)
!
!.....Call the related subroutine for the calculation of Cd
! 
      ncell_fluid=nag*nt
      CALL int_drag
      Cdl(1:ncell_fluid)=cell%vfgl(1:ncell_fluid)
      Cdd(1:ncell_fluid)=cell%vfgd(1:ncell_fluid)
!
!.....Write output files
!
      DO i=1,nag
         WRITE(800,100) cell%alphag(i),(H_il(i+nag*(j-1)),j=1,nt)
         WRITE(900,100) cell%alphag(i),(H_ig(i+nag*(j-1)),j=1,nt)
      ENDDO
!
      ii=0
      DO i=1,nt
         ii=1+nag*(i-1)
         j=nag*i
         WRITE(801,100) cell%tl(ii),H_il(ii),H_il(j)
         WRITE(901,100) cell%tg(ii),H_ig(ii),H_ig(j)
      ENDDO
!
      ii=0
      WRITE(802,*) 'zone i=',nag,',j=',nt
      WRITE(902,*) 'zone i=',nag,',j=',nt
      DO j=1,nt
         DO i=1,nag
            ii=ii+1
            WRITE(802,100) cell%tl(ii),cell%alphag(ii),H_il(ii)
            WRITE(902,100) cell%tg(ii),cell%alphag(ii),H_ig(ii)
         ENDDO
      ENDDO
!
      ii=0
      WRITE(803,*) 'zone i=',nag,',j=',nt
      WRITE(903,*) 'zone i=',nag,',j=',nt
      DO j=1,nt
         DO i=1,nag
            ii=ii+1
            WRITE(803,100) gamma_void(ii),cell%alphag(ii),H_il(ii)
            WRITE(903,100) gamma_void(ii),cell%alphag(ii),H_ig(ii)
         ENDDO
      ENDDO
!
      ii=0
      WRITE(804,*) 'zone i=',nag,',j=',nt
      WRITE(904,*) 'zone i=',nag,',j=',nt
      DO j=1,nt
         DO i=1,nag
            ii=ii+1
            WRITE(804,100) cell%tl(ii),cell%alphag(ii),dHldag(ii)
            WRITE(904,100) cell%tg(ii),cell%alphag(ii),dHgdag(ii)
         ENDDO
      ENDDO
!
      ii=0
      WRITE(805,*) 'zone i=',nag,',j=',nt
      WRITE(905,*) 'zone i=',nag,',j=',nt
      DO j=1,nt
         DO i=1,nag
            ii=ii+1
            WRITE(805,100) cell%tl(ii),cell%alphag(ii),dHldtl(ii)
            WRITE(905,100) cell%tg(ii),cell%alphag(ii),dHgdtg(ii)
         ENDDO
      ENDDO
!
      ii=0
      WRITE(810,*) 'zone i=',nag,',j=',nt
      DO j=1,nt
         DO i=1,nag
            ii=ii+1
            WRITE(810,100) gamma_void(ii),cell%alphag(ii),Cdl(ii)
         ENDDO
      ENDDO
!
  100 FORMAT(30(e16.5,1x))
!
      PAUSE 'Check the Cd, Hik surface..'
!
!.....Terminate the program after writing Hik
!
      STOP
      END SUBROUTINE udfn_check_Hik
