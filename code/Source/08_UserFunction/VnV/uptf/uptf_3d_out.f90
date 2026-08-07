!
      SUBROUTINE uptf_3D_out    
!
!     This routine controls user/problem specific output variables & format
!
      USE VOL_DATA        , ONLY: cell
      USE Zzone           , ONLY: ncell_fluid,ncell_fluid_all
      USE Zb_condition    , ONLY: alphab_liq,alphab_gas,rhob_liq,rhob_gas
      USE Zconst2         , ONLY: iprn
      USE Zcoord3         , ONLY: vol
      USE Zcore           , ONLY: np,myrank  
!     USE Zflux           , ONLY: fluxvol_g,fluxvol_l      
      USE Zvec_major      , ONLY: flux_l_nf,flux_g_nf
      USE Znum_cell       , ONLY: istart_nf,istart_nbcon_nf
      USE Ztimecon        , ONLY: time,itim
      USE Zuptf           , ONLY: cwl
      USE Zvec_index    , ONLY: left_nf,nbcon_nf
!      
      IMPLICIT NONE
!
      INTEGER :: i,k,ii,na
      INTEGER :: nf_number,istart,len,istart2,i1,i2
!
!.....Local variables
      REAL(8) :: ECC_deliv,bypass                           !for UPTF 
      REAL(8) :: tmp(4)
      REAL(8) :: level      
      REAL(8) :: in1,in2,in3,in4,in5,in6,out1,out2  !for UPTF
      REAL(8),SAVE :: CALL_time
      REAL(8),SAVE :: intv_call=0.1d0
!.....Local allocatable arrays
      REAL(8),DIMENSION(:),ALLOCATABLE :: vdt_dat1
!
      LOGICAL, SAVE::initial=.true.
!
      na=ncell_fluid_all
      IF(initial) THEN
         initial=.false.
         CALL_time=time    !time for data saving for averaging calc.
         IF(myrank.eq.0) THEN
            OPEN(130,file='UPTF_balance.dat')  
            OPEN(888,file='UPTF_collapsed_water_level.dat')
         ENDIF
      ENDIF
         
      IF(time.ge.CALL_time)THEN
         CALL_time=CALL_time+intv_call      
         in1=0.0d0
         in2=0.0d0
         in3=0.0d0
         in4=0.0d0
         in5=0.0d0
!            
         nf_number=2
         istart=istart_nf(1,nf_number)
         istart2=istart_nbcon_nf(nf_number)
         len   =istart_nf(2,nf_number)
         DO i=1,len
            i1=istart+i
            i2=istart2+i
            ii=left_nf(i1)
            k=nbcon_nf(i2)
            IF(k.eq.1) in1=in1+alphab_liq(k)*rhob_liq(k)*flux_l_nf(i1)+alphab_gas(k)*rhob_gas(k)*flux_g_nf(i1)
            IF(k.eq.2) in2=in2+alphab_liq(k)*rhob_liq(k)*flux_l_nf(i1)+alphab_gas(k)*rhob_gas(k)*flux_g_nf(i1)
            IF(k.eq.3) in3=in3+alphab_liq(k)*rhob_liq(k)*flux_l_nf(i1)+alphab_gas(k)*rhob_gas(k)*flux_g_nf(i1)
            IF(k.eq.4) in4=in4+alphab_liq(k)*rhob_liq(k)*flux_l_nf(i1)+alphab_gas(k)*rhob_gas(k)*flux_g_nf(i1)
            IF(k.eq.5) in5=in5+alphab_liq(k)*rhob_liq(k)*flux_l_nf(i1)+alphab_gas(k)*rhob_gas(k)*flux_g_nf(i1)
            IF(k.eq.5) in6=cell%p(i)
         ENDDO
!            
         out1=0.d0
         out2=0.d0
         nf_number=3
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
         DO i=1,len
            i1=istart+i
            ii=left_nf(i1)
            out1=cell%alphal(ii)*cell%rhol(ii)*flux_l_nf(i1)+cell%alphag(ii)*cell%rhog(ii)*flux_g_nf(i1)
            out2=cell%alphag(ii)
         ENDDO
!
         IF(np.gt.1) THEN
            tmp(1)=ABS(in1)
            tmp(2)=ABS(in2)
            tmp(3)=ABS(in3)
            tmp(4)=ABS(in4)
            CALL reducei_max_r(tmp,4)
            in1=tmp(1)
            in2=tmp(2)
            in3=tmp(3)
            in4=tmp(4)
!
            tmp(1)=ABS(in5)
            tmp(2)=ABS(out1)
            tmp(3)=ABS(out2)
            CALL reducei_r(tmp,3)
            in5=tmp(1)
            out1=tmp(2)
            out2=tmp(3)
         ELSE
            in1=ABS(in1)
            in2=ABS(in2)
            in3=ABS(in3)
            in4=ABS(in4)
            in5=ABS(in5)
            out1=DABS(out1)
            out2=DABS(out2)
         ENDIF
!
         IF(myrank.eq.0) THEN
            bypass=out1-in5
            ECC_deliv=in1+in2+in3-bypass
            WRITE(130,1130)time,ECC_deliv,bypass,in1,in2,in3,in4,in5,out1,in6,out2
         ENDIF
!             
      ENDIF
!
!.....collapsed water lever        
!
!     
      tmp(1)=0.d0
      tmp(2)=0.d0
      DO i=1,ncell_fluid
         tmp(1) = tmp(1) + cell%alphal(i)*cell%rhol(i)*vol(i)
         tmp(2) = tmp(2) + cell%rhol(i)*vol(i)
      ENDDO
      IF(np.gt.1) CALL allreducei_r(tmp,2)
      level = tmp(1)/tmp(2)*10.2d0 
      cwl=level
!            
      IF(MOD(itim,iprn).eq.0)THEN
         IF(myrank.eq.0) WRITE(888,250) time, level  
      ENDIF   
!      
!.....Compare np1 vs np4         
      IF(myrank.eq.0) THEN
         ALLOCATE(vdt_dat1(na))
      ELSE
         ALLOCATE(vdt_dat1(1))
      ENDIF
      CALL gatherv_r(cell%alphal,ncell_fluid,vdt_dat1,na,0)
      IF(np.eq.1) THEN
         OPEN(431, file='VFTX_ref-compare_np1.dat')
         DO i=1,na
            WRITE(431,311) i, vdt_dat1(i)
         ENDDO                  
         CLOSE(431)
      ELSE
         IF(myrank.eq.0) THEN
            OPEN(431, file='VFTX_ref-compare_np4.dat')
            DO i=1,na
               WRITE(431,311) i, vdt_dat1(i)
            ENDDO               
            CLOSE(431)
         ENDIF   
      ENDIF   
      DEALLOCATE(vdt_dat1)  
!
311   FORMAT(i5,2x,e14.7)     
250   FORMAT(4x,4e15.7)           
1130  FORMAT(1x,30f20.7)    
!
      END SUBROUTINE  uptf_3D_out    
