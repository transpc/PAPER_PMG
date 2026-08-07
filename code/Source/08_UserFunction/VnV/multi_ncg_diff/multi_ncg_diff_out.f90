!
      SUBROUTINE multi_ncg_diff_out
!
      USE Zinterface
      USE VOL_DATA     , ONLY: cell
      USE Zzone        , ONLY: ncell_fluid,ncell_fluid_all
      USE Zcore        , ONLY: myrank
      USE Zvec_param   , ONLY: nf_out
      USE Znum_cell    , ONLY: istart_nf,nf_number_nb,lens,       &
                               istart_nfs,nf_number_id,istart_nf
      USE Ztimecon     , ONLY: time
      USE Zconst2      , ONLY: dt
      USE Zbc_index    , ONLY: npb
      USE Zb_condition , ONLY: eb_gas,rhob_gas,tb_gas,vb_gas
      USE Zcoord3      , ONLY: volp
      USE Zncg         , ONLY: qn_cell,ra_nvin
      USE Zvec_major   , ONLY: flux_g_nf
!
      IMPLICIT NONE
!
      INTEGER,PARAMETER :: nout=10
!.....Local variables
      INTEGER i,na
      INTEGER :: nv,nf_number,len,istart0,istart,i0,i1
      INTEGER,SAVE :: iout=0
      LOGICAL,SAVE:: initial=.true.
      LOGICAL,SAVE:: init=.true.
      REAL(8) :: tot_eg,tot_hg,eg_err,hg_err
      REAL(8) :: tot_mass,mass_err,ar
      REAL(8),SAVE :: tot_eg0,tot_hg0,eg_in,eg_out,hg_in,hg_out
      REAL(8),SAVE :: tot_mass0,mass_in,mass_out
!.....Local arrays
      REAL(8),DIMENSION(ncell_fluid) :: flux_g
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: dat
!.....Local vector arrays
      REAL(8),DIMENSION(nf_out) :: flux_g_out
!
      IF(initial.and.myrank.eq.0)THEN
         OPEN(990,file='VFT13_multi_ncg_ref.dat')
         initial=.false.
      ENDIF
!
      tot_eg  =0.d0
      tot_hg  =0.d0
      tot_mass=0.d0
      DO i=1,ncell_fluid
         IF(npb(i).eq.0)THEN
            tot_eg  =tot_eg  +cell%rhog(i)*cell%eg(i)*volp(i)
            tot_hg  =tot_hg  +cell%rhog(i)*cell%ha(i)*volp(i)
            tot_mass=tot_mass+cell%rhog(i)*volp(i)
         ENDIF
      ENDDO
!
      IF(init)THEN
         tot_eg0  =tot_eg
         tot_hg0  =tot_hg
         tot_mass0=tot_mass
         eg_in   =0.d0
         eg_out  =0.d0
         hg_in   =0.d0
         hg_out  =0.d0
         mass_in =0.d0
         mass_out=0.d0
         init=.false.
      ENDIF
!
      ar=1.d-2
!      ar=0.0016d0
!
      eg_in  =eg_in  +vb_gas(1,3)*ar* eb_gas(1)                      *rhob_gas(1)*dt
      hg_in  =hg_in  +vb_gas(1,3)*ar*(eb_gas(1)+ra_nvin(1)*tb_gas(1))*rhob_gas(1)*dt
      mass_in=mass_in+vb_gas(1,3)*ar                                 *rhob_gas(1)*dt
!
!.....Build summation info for out
!
      nf_number_nb=0
      nf_number_id(0)=3
      istart_nfs(0)=0
      lens         =istart_nfs(0) +nf_out
!
      nv=0
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      IF(len.gt.0) THEN
         DO i=1,len
            i0=istart0+i
            i1=istart+i
            flux_g_out(i0)=flux_g_nf(i1)
         ENDDO
!
         CALL sum_nf(0,1,               &
                     flux_g_out,flux_g)
!
         DO i=1,ncell_fluid
            eg_out  =eg_out  +flux_g(i)*cell%eg(i)*cell%rhog(i)
            hg_out  =hg_out  +flux_g(i)*cell%ha(i)*cell%rhog(i)
            mass_out=mass_out+flux_g(i)*cell%rhog(i)
         ENDDO
      ENDIF
      eg_out  =eg_out*dt
      hg_out  =hg_out*dt
      mass_out=mass_out*dt
!
      eg_err  =(tot_eg0+eg_in-eg_out-tot_eg)        /tot_eg0
      hg_err  =(tot_hg0+hg_in-hg_out-tot_hg)        /tot_hg0
      mass_err=(tot_mass0+mass_in-mass_out-tot_mass)/tot_mass0
!
      iout=iout+1
      IF(iout.ge.nout)THEN
         na=ncell_fluid_all
         IF(myrank.eq.0) THEN
            ALLOCATE(dat(na,3))
         ELSE
            ALLOCATE(dat(1,3))
         ENDIF
         CALL gatherv_r(cell%tg     ,ncell_fluid,dat(1,1),na,0)
         CALL gatherv_r(qn_cell(1,1),ncell_fluid,dat(1,2),na,0)
         CALL gatherv_r(qn_cell(1,2),ncell_fluid,dat(1,3),na,0)
!         IF(myrank.eq.0) WRITE(990,10) time,dat(14,3),dat(14,2),dat(14,1),eg_err,hg_err,mass_err,tot_eg,tot_hg
         IF(myrank.eq.0) THEN
            WRITE(990,10) time,dat(14,3),dat(14,2),dat(14,1)
         ENDIF
         iout=0
         DEALLOCATE(dat)
      ENDIF
!
   10 FORMAT(20(4e15.7,1x))
!
      END SUBROUTINE multi_ncg_diff_out
