!
      SUBROUTINE wall_condensation_steam
!
!     This subroutine calculate wall condensation for only steam.
!     Liquid properties are defined by face temperature.
!
      USE Zinterface
      USE VOL_DATA        , ONLY: cell
      USE Wall_DATA       , ONLY: face
      USE Zvec_param      , ONLY: nf_ctw,nf_fsw
      USE Zvec_index      , ONLY: left_nf
      USE Znum_cell       , ONLY: istart_nf,nf_number_nb,lens,nf_number_id,istart_nfs
      USE Zcoord3         , ONLY: volr
      USE Zqvol           , ONLY: gamma_wall
      USE Zmodel          , ONLY: qconden
      USE Zvec_index_solid, ONLY: qliq_fsw,qgas_fsw,vfilm_fsw,dfilm_fsw,vfilm_ctw,dfilm_ctw
      USE Zvec_geo        , ONLY: saa_nf
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i,ii
      INTEGER :: nv,nf_number,len,istart0,istart,i0,i1
      REAL(8) :: hcond
!.....Local vector arrays
      REAL(8),DIMENSION(nf_fsw) :: hcond_fsw
      REAL(8),DIMENSION(nf_ctw) :: hcond_ctw,twall_ctw
      REAL(8),DIMENSION(nf_fsw+nf_ctw) :: gamma_wall_nf,qconden_nf
!
!.....Build summation info for fsw,ctw
!
      nf_number_nb=1
      nf_number_id(0)=5
      nf_number_id(1)=6
      istart_nfs(0)=0
      istart_nfs(1)=istart_nfs(0)+nf_fsw
      lens         =istart_nfs(1)+nf_ctw
!
!.....Calculate condensation heat transfer coefficient along the face:
!.....constant temperature face or fluid-solid interface.
!
!.....Cells fsw
!
      nv=0
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      IF(nf_fsw.gt.0) CALL wall_condensation_hcond(hcond_fsw,dfilm_fsw,vfilm_fsw,nf_number,nf_fsw)
      DO i=1,len  
         i0=istart0+i
         i1=istart+i
         ii=left_nf(i1)
         IF(face%twall_partition(ii).lt.cell%ts(ii))THEN
            hcond=hcond_fsw(i)
            qconden_nf(i0)=-hcond*(cell%ts(ii)-face%twall_partition(ii))
            gamma_wall_nf(i0)=qconden_nf(i0)*saa_nf(i1)*volr(ii)/(cell%hg(ii)-cell%hlsat(ii))            
            qliq_fsw(i)=0.0d0
            qgas_fsw(i)=0.0d0
         ENDIF
      ENDDO
!
!.....Cells ctw
!
      CALL udfn_tw_profile(twall_ctw)
      nv=1
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      DO i=1,len
         i1=istart+i
         ii=left_nf(i1)
         face%twall_partition(ii)=twall_ctw(i)
      ENDDO
      IF(nf_ctw.gt.0) CALL wall_condensation_hcond(hcond_ctw,dfilm_ctw,vfilm_ctw,nf_number,nf_ctw)
      DO i=1,len  
         i0=istart0+i
         i1=istart+i
         ii=left_nf(i1)
         IF(face%twall_partition(ii).lt.cell%ts(ii))THEN
            hcond=hcond_ctw(i)
            qconden_nf(i0)=-hcond*(cell%ts(ii)-face%twall_partition(ii))
            gamma_wall_nf(i0)=qconden_nf(i0)*saa_nf(i1)*volr(ii)/(cell%hg(ii)-cell%hlsat(ii))
         ENDIF
      ENDDO
!
      CALL sum_nf(0,-1,                     &
                  gamma_wall_nf,gamma_wall, &
                  qconden_nf   ,qconden)
!
      END SUBROUTINE wall_condensation_steam  
!
      SUBROUTINE wall_condensation_hcond(hcond,dfilm,vfilm,nf_number,nf_len)
!
      USE VOL_DATA        , ONLY: cell
      USE Wall_DATA       , ONLY: face
      USE Zvec_index      , ONLY: left_nf
      USE Znum_cell       , ONLY: istart_nf
      USE STM_TBL_cupid   , ONLY: st_tbl,    &
                                  nt,ndxstd, &
                                  pcrit
      USE Zmodel          , ONLY: wVertical
!
      IMPLICIT NONE
!
!.....Input
      REAL(8) dfilm(nf_len),vfilm(nf_len)
      INTEGER::nf_number,nf_len
!.....Output
      REAL(8) hcond(nf_len)
!.....Local variables
      INTEGER i,ii
      INTEGER istart,len,i1
      LOGICAL err
      REAL(8) Tliq(1),Rhof(1),Cpf,Muf(1),Kf(1),vis_rho
      REAL(8) htdiam,Vliq,refilm,filmt,hfg,twsub
      REAL(8) quax,rey,pr,hdb,hf,z,ftr
      REAL(8) hshah
      REAL(8) t1
!.....Local arrays
      REAL(8) s(36)
!
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      DO i=1,len  
         i1=istart+i
         ii=left_nf(i1)
!         
! Three 'face%twallpartitions' exist in user_def_inp, scalar_eneryg_diffusion and in scalar_matrix_solid,
! and you must see which value is used finally.
! Twall case uses that of user_def_inp, and
! F-S case uses that of scalar_energy_diffusion or scalar_matrix_solid.
!         
         IF(face%twall_partition(ii).lt.cell%ts(ii))THEN
            Tliq(1)=face%twall_partition(ii)
            s(1)=Tliq(1)
            s(9)=0.0d0
            CALL sth2x1_cupid(st_tbl(ndxstd),st_tbl(ndxstd+nt),s,err)
            Rhof(1)=1.d0/s(11)
            Cpf=s(21)
            CALL viscos_lw_cupid(Tliq,Rhof,cell%ts(ii),Muf,'liq',1)
            CALL cond_lw_cupid(Tliq,Rhof,Kf,1)
!
            htdiam=dfilm(i) ! <-- input
            Vliq=vfilm(i)   ! <-- input
!            htdiam=0.0005d0       
!            Vliq=SQRT(vl_n(ii,1)**2+vl_n(ii,2)**2+vl_n(ii,3)**2)
!
!...........Nusselt : laminar film, vertical tube/wall
!
            IF(wVertical(ii).eq.1)THEN
               vis_rho=Muf(1)/Rhof(1)
               refilm=Rhof(1)*Vliq*htdiam/Muf(1)   ! (MARS ; gliqa??)
               filmt=0.9086d0*(refilm*vis_rho*vis_rho/9.81d0)**0.33333d0
               filmt=MAX(1.0d-05,filmt)      ! limit film thickness : > 10 microns
               hcond(i)=Kf(1)/filmt
!
!...........Chato : laminar film, horizontal tube
!
            ELSEIF(wVertical(ii).eq.-1)THEN
               hfg=cell%hg(ii)-cell%hlsat(ii)
               twsub=cell%ts(ii)-face%twall_partition(ii)
               t1=sqrt(Rhof(1)*Rhof(1)*9.81d0*hfg*Kf(1)*Kf(1)*Kf(1)/(htdiam*Muf(1)*twsub))
               hcond(i)=0.296d0*t1*sqrt(t1)
            ENDIF
!
!...........shah : turbulence film, vertical/horizontal
!
!...........pcrit is constant read in stread.f90 why change?
!           pcrit=22.4d6
!           quax=cell%alphag(ii)            ! (MARS ; v_da(iv)%Voidg -> v_da(iv)%Quale -> quax)
            quax=max(min(cell%alphag(ii),1.0d0),1.0d-9)
            rey=Rhof(1)*Vliq*htdiam/Muf(1)         ! (MARS ; gabs??)
            pr=Muf(1)*Cpf/Kf(1)
            hdb=0.023d0*Kf(1)*rey**0.8d0*pr**0.4d0/htdiam
            hf=hdb*(1.0d0-quax)**0.8d0       !quax=1.0 always
            z=(cell%p(ii)/pcrit)**0.4d0*(1.0d0/quax-1.0d0)**0.8d0
            ftr=1.0d0
            IF(z.ne.0.0d0) ftr=1.0d0+3.8d0/z**0.95d0 
            hshah=hf*ftr
            hcond(i)=MAX(hcond(i),hshah)     !hshah=0.0d0 due to quax=1.0d0 always
         ENDIF
      ENDDO
!
      END SUBROUTINE wall_condensation_hcond
