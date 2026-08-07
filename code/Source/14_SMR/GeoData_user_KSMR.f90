      subroutine GeoData_user_KSMR(nn)
!
!      geometry info of KSMR 
!      nn=ncell --> global cell number (total)
!    
    
      USE Zcore       , ONLY: myrank,np
      USE Zrv_subchan  , ONLY: subchannel_type_tmp
      USE Zrv_ncell   , ONLY: asm_nx,asm_ny,asm_nz,asm_ni,asm_ni2,chn_nx,chn_ny
      USE Zrv_hts_2d  , ONLY: dz_fuel00
      USE Zporous     , ONLY: chn_type_tmp
   
      INTEGER i,ii,size,icell,icnt,nch_FA,num_row,iFA,nch,ilev,num_FA,id,istart,ichx,ichy
      INTEGER nz0,nn
      INTEGER, ALLOCATABLE :: rowFA(:)
      
      REAL(8) :: dheight
!
!User input    
      nch=17          !16x16 rod bundle/FA = 17x17 channels/FA 
      nch_FA=nch*nch   !number of channels in each FA
      num_FA=69       !number of FA
      num_row=9       !number of FA rows
      nz0=10          !number divisions of Core in axial direction 
      dheight=0.25d0   !cell height 
!         
      ALLOCATE(rowFA(num_row))
!    
      ALLOCATE(asm_nx(nn),asm_ny(nn),asm_nz(nn),asm_ni(nn),asm_ni2(nn)) 
      ALLOCATE(chn_nx(nn),chn_ny(nn))
      ALLOCATE(subchannel_type_tmp(nn))
      ALLOCATE(chn_type_tmp(nn))
!
      DO i=1,nn  !nn=ncell
         asm_nx(i)=0
         asm_ny(i)=0
         asm_nz(i)=0
         asm_ni(i)=0
         asm_ni2(i)=0
         chn_nx(i)=0
         chn_ny(i)=0
         subchannel_type_tmp(i)=0
         chn_type_tmp(i)=0
      ENDDO    
!
      IF(myrank.eq.0) then
         OPEN(unit=1551,file='core.dat')     !to read
         OPEN(unit=1999,file='core_mod.dat')  !to check  
         READ(1551,*) size
!       
         icnt=0
         iFA=0
         ichx=0
         ichy=1
         ilev=0
         DO ii=1,size
            icnt=icnt+1 !nrod
            READ(1551,*) icell
!
!User input          
!FA geometry information (Sum of FAs= end-start+1 of each FA row)
            rowFA(1)=7-3+1          !number of FAs of the 1st row (at present 5 FAs in the 1st row) 
            rowFA(2)=rowFA(1)+8-2+1
            rowFA(3)=rowFA(2)+9-1+1
            rowFA(4)=rowFA(3)+9-1+1
            rowFA(5)=rowFA(4)+9-1+1
            rowFA(6)=rowFA(5)+9-1+1
            rowFA(7)=rowFA(6)+9-1+1
            rowFA(8)=rowFA(7)+8-2+1
            rowFA(9)=rowFA(8)+7-3+1
!
!User input          
!FA row id and starting FA number          
            IF(icnt.le.nch_FA*rowFA(1)) THEN
               id=1
               istart=3
            ELSEIF(icnt.le.nch_FA*rowFA(2)) THEN    
               id=2
               istart=2
            ELSEIF(icnt.le.nch_FA*rowFA(3)) THEN    
               id=3
               istart=1
            ELSEIF(icnt.le.nch_FA*rowFA(4)) THEN    
               id=4
               istart=1
            ELSEIF(icnt.le.nch_FA*rowFA(5)) THEN    
               id=5
               istart=1
            ELSEIF(icnt.le.nch_FA*rowFA(6)) THEN    
               id=6
               istart=1
            ELSEIF(icnt.le.nch_FA*rowFA(7)) THEN    
               id=7
               istart=1
            ELSEIF(icnt.le.nch_FA*rowFA(8)) THEN    
               id=8
               istart=2
            ELSE
               id=9
               istart=3
            ENDIF    
!          
            IF(icnt.le.nch_FA*rowFA(id)) THEN
!
!Define Assemby coordinate              
               asm_nx(icell)=istart+iFA
               asm_ny(icell)=id
               asm_nz(icell)=1+ilev
               IF(MOD(icnt,nch_FA).eq.0) iFA=iFA+1
               IF(MOD(icnt,nch_FA*rowFA(id)).eq.0) iFA=0
!              
!Define channel coordinate              
               ichx=ichx+1              
               chn_nx(icell)=ichx
               chn_ny(icell)=ichy
               IF(ichx.eq.nch) THEN
                  ichx=0
                  ichy=ichy+1
               ENDIF
               IF(MOD(icnt,nch_FA).eq.0) THEN
                  ichx=0
                  ichy=1
               ENDIF    
!
!User input              
!Define subchannel type              
!                center              
               subchannel_type_tmp(icell)=1
               chn_type_tmp(icell)=1
!              corner              
               IF( (chn_nx(icell).eq.1.and.chn_ny(icell).eq.1)     .or. &
                   (chn_nx(icell).eq.1.and.chn_ny(icell).eq.nch)   .or. &
                   (chn_nx(icell).eq.nch.and.chn_ny(icell).eq.1)   .or. &
                   (chn_nx(icell).eq.nch.and.chn_ny(icell).eq.nch)       ) THEN
                   subchannel_type_tmp(icell)=3
                   chn_type_tmp(icell)=3
               ENDIF   
!              side
               IF( (chn_nx(icell).gt.1.and.chn_nx(icell).lt.nch.and.chn_ny(icell).eq.1)     .or. &
                   (chn_nx(icell).gt.1.and.chn_nx(icell).lt.nch.and.chn_ny(icell).eq.nch)   .or. &
                   (chn_ny(icell).gt.1.and.chn_ny(icell).lt.nch.and.chn_nx(icell).eq.1)     .or. &
                   (chn_ny(icell).gt.1.and.chn_ny(icell).lt.nch.and.chn_nx(icell).eq.nch)           ) THEN
                   subchannel_type_tmp(icell)=2
                   chn_type_tmp(icell)=2
               ENDIF               
!             Test
              !subchannel_type_tmp(icell)=0 
               
           ENDIF    
           IF(MOD(icnt,nch_FA*num_FA).eq.0) THEN
              ilev=ilev+1
              icnt=0
           ENDIF    
!          
           WRITE(1999,199) icell,asm_nx(icell),asm_ny(icell),asm_nz(icell),chn_nx(icell),chn_ny(icell),subchannel_type_tmp(icell)
         ENDDO
         
         CLOSE(1551)
         CLOSE(1999)         
      ENDIF  
199   FORMAT(100(i9,3x))                 
!
      IF(np.gt.1) THEN
         CALL broadcast_i(asm_nx,nn)
         CALL broadcast_i(asm_ny,nn)
         CALL broadcast_i(asm_nz,nn)
         CALL broadcast_i(asm_ni,nn)
         CALL broadcast_i(asm_ni2,nn)
         CALL broadcast_i(chn_nx,nn)
         CALL broadcast_i(chn_ny,nn)
         CALL broadcast_i(subchannel_type_tmp,nn)
         CALL broadcast_i(chn_type_tmp,nn)
      ENDIF   
!
!.....Axial 
!
      ALLOCATE(dz_fuel00(nz0))
      IF(myrank.eq.0) THEN
         !DO j=1,nz0
         !   READ(1550,*)dz_fuel00(j)
         !ENDDO
         dz_fuel00(1:nz0)=0.24d0
      ENDIF
      IF(np.gt.1) CALL broadcast_i1(nz0)
      IF(np.gt.1) CALL broadcast_r(dz_fuel00,nz0)
!
    RETURN
    END SUBROUTINE
    