!
      SUBROUTINE rocom_outlet_bc_user
!
!     Define an outlet boundary conditoin to the ROCOM problem
!
      USE VOL_DATA        , ONLY: cell
      USE Zparam          , ONLY: nin_max,nb_max
      USE Zvec_param      , ONLY: nf_inl,nf_out
      USE Ztimecon        , ONLY: time
      USE Znum_cell       , ONLY: istart_nf,istart_nbcon_nf
      USE Zvec_index      , ONLY: left_nf,nbcon_nf
      USE Zb_condition    , ONLY: pbnd
      USE Zpress          , ONLY: p
      USE Zrocom_specific , ONLY: mdot_bc,out_mdot,in_mdot
      USE Zbc_index       , ONLY: npin
      USE Zvec_major      , ONLY: flux_l_nf,flux_g_nf,flux_d_nf
!
      IMPLICIT NONE      
!      
!.....Local variables
      INTEGER :: i,j,k
      INTEGER :: nf_number,istart,isize,istart2,i1,i2
      INTEGER :: ii,nmax_pin
      LOGICAL,SAVE :: INITIAL=.TRUE.
!.....Local arrays
      INTEGER :: ipbnd_x(npin)
      REAL(8) :: pbnd_x(npin)
      REAL(8) :: mdot_bc_inl(nf_inl),mdot_bc_out(nf_out)
      INTEGER,DIMENSION(:),SAVE,ALLOCATABLE :: obr_loop,obr_loopmax,optpupDOwn
      REAL(8),DIMENSION(:),SAVE,ALLOCATABLE :: pupdown
!
      IF(1 .and. time.gt.0.0d0)THEN      
!
!........Cells inl
!
         nf_number=2
         istart=istart_nf(1,nf_number)
         isize =istart_nf(2,nf_number)
         DO i=1,isize
            i1=istart+i
            ii=left_nf(i1)
            mdot_bc_inl(i)= cell%alphal(ii)*cell%rhol(ii)*flux_l_nf(i1) &
                           +cell%alphag(ii)*cell%rhog(ii)*flux_g_nf(i1) &
                           +cell%alphad(ii)*cell%rhod(ii)*flux_d_nf(i1)
         ENDDO
!
!........Cells out
!
         nf_number=3
         istart=istart_nf(1,nf_number)
         isize =istart_nf(2,nf_number)
         DO i=1,isize
            i1=istart+i
            ii=left_nf(i1)
            mdot_bc_out(i)= cell%alphal(ii)*cell%rhol(ii)*flux_l_nf(i1) &
                           +cell%alphag(ii)*cell%rhog(ii)*flux_g_nf(i1) &
                           +cell%alphad(ii)*cell%rhod(ii)*flux_d_nf(i1)
         ENDDO
!     
         CALL sum_nf23_nbcon(mdot_bc_inl,mdot_bc_out,mdot_bc)
!
!........sum_nf23_nbcon takes care of minus sign
!
!        DO i=1, nin_max
!           mdot_bc(i)=-mdot_bc(i)
!        ENDDO
         DO i=1,nin_max   
            in_mdot(i)=mdot_bc(i)
         ENDDO
         DO i=nin_max+1,nb_max
            out_mdot(i-nin_max)=mdot_bc(i)
         ENDDO
!     
!........At the first CALL of this SUBROUTINE      
!
         IF(INITIAL)THEN
            INITIAL=.FALSE.
            ALLOCATE(pupdown(nb_max))
            ALLOCATE(obr_loop(nb_max),obr_loopmax(nb_max),optpupDOwn(nb_max))
            pupDOwn(:)      =1.d0 
            pupDOwn(1:4)    =5. !pa
            optpupDOwn(1:4) =0
            obr_loop(1:4)   =0
            obr_loopmax(1:4)=1
         ENDIF 
!        
!........Adjust pressure at 8th outlet(4th npb)  
!
         nmax_pin=0
         DO i=1,npin
            IF(obr_loop(i).gt.obr_loopmax(i))THEN
               obr_loop(i)=0
               IF(out_mdot(i).gt.in_mdot(i)*1.2)THEN
                  IF(pupdown(i).eq.1) pupdown(i)=pupdown(i)*1.2d0
                  IF(pupdown(i).eq.-1) pupdown(i)=pupdown(i)*0.9d0
                  optpupdown(i)=1
               ELSEIF(out_mdot(i).gt.in_mdot(i)*1.1)THEN
                  IF(pupdown(i).eq.1) pupdown(i)=pupdown(i)*1.1d0
                  IF(pupdown(i).eq.-1) pupdown(i)=pupdown(i)*0.95d0
                  optpupdown(i)=1
               ELSEIF(out_mdot(i).gt.in_mdot(i)*1.05)THEN
                  IF(pupdown(i).eq.1) pupdown(i)=pupdown(i)*1.05d0
                  IF(pupdown(i).eq.-1) pupdown(i)=pupdown(i)*0.975d0
                  optpupdown(i)=1
               ELSEIF(out_mdot(i).gt.in_mdot(i)*1.01)THEN
                  optpupdown(i)=1
               ELSEIF(out_mdot(i).gt.in_mdot(i)*1.001)THEN
                  optpupdown(i)=1   
               ELSEIF(out_mdot(i).lt.in_mdot(i)*0.8)THEN
                  IF(pupdown(i).eq.-1) pupdown(i)=pupdown(i)*1.2d0
                  IF(pupdown(i).eq.1) pupdown(i)=pupdown(i)*0.9d0
                  optpupdown(i)=-1
               ELSEIF(out_mdot(i).lt.in_mdot(i)*0.9)THEN
                  IF(pupdown(i).eq.-1) pupdown(i)=pupdown(i)*1.1d0
                  IF(pupdown(i).eq.1) pupdown(i)=pupdown(i)*0.95d0
                  optpupdown(i)=-1
               ELSEIF(out_mdot(i).lt.in_mdot(i)*0.95)THEN
                  IF(pupdown(i).eq.-1) pupdown(i)=pupdown(i)*1.05d0
                  IF(pupdown(i).eq.1) pupdown(i)=pupdown(i)*0.975d0
                  optpupdown(i)=-1
               ELSEIF(out_mdot(i).lt.in_mdot(i)*0.99)THEN
                  optpupdown(i)=-1
               ELSEIF(out_mdot(i).lt.in_mdot(i)*0.999)THEN
                  optpupdown(i)=-1
               ELSE 
                  optpupdown(i)=0            
               ENDIF
               pbnd(i)=pbnd(i)+optpupdown(i)*pupdown(i)
!
!..............Limit pressure of the pressure boundary      
!
               IF(pbnd(i).gt.38.d5) pbnd(i)=38.d5
               IF(pbnd(i).lt.37.990d5) pbnd(i)=37.990d5
!! save the i needed to update p(icell)
               nmax_pin=nmax_pin+1
               ipbnd_x(nmax_pin)=i
               pbnd_x(nmax_pin)=pbnd(i)
!!
!               
!..............Set the pressure of the pressure cell
!
! not efficient to scan nbcon array for each found i 
!              DO icell=1,ncell_fluid
!                 j0=i_neigh(icell)-1
!                 DO jth=i_neigh(icell),i_neigh(icell+1)-1
!                    IF(nbcon(jth-j0,icell).eq.nin_max+i) p(icell)=pbnd(i)
!                 ENDDO
!              ENDDO         
!               
            ENDIF
            obr_loop(i)=obr_loop(i)+1 
         ENDDO
!!
         nf_number=3
         istart=istart_nf(1,nf_number)
         istart2=istart_nbcon_nf(nf_number)
         isize =istart_nf(2,nf_number)
         DO i=1,isize
            i1=istart+i
            i2=istart2+i
            ii=left_nf(i1)
            DO k=1,nmax_pin
               j=ipbnd_x(k)
               IF(nbcon_nf(i2).eq.nin_max+j) p(ii)=pbnd_x(k)
            ENDDO
         ENDDO
      ENDIF
!
      END SUBROUTINE rocom_outlet_bc_user
