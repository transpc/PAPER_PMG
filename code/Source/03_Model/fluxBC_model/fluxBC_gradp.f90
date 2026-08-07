      SUBROUTINE fluxBC_gradp(s,fie_nf)

      USE Zmpi         , ONLY: ncell_fp
      USE Zzone        , ONLY: ncell_fluid
      USE Zparam       , ONLY: ndim    
      USE Zvec_param   , ONLY: nf_totk
      USE Zvec_index   , ONLY: left_nf,right_non      
      USE Znum_cell    , ONLY: istart_nf,istart_nfs,nf_number_id,right_nb_k
      USE Zvec_geo     , ONLY: sv_nf
      USE Zrv_choke    , ONLY: num_throatface,n_face_throat,choke,nonk_throat
      USE Zvalve       , ONLY: valve_closed,num_valveloc,num_valveface,n_face_valve,mapping_valve,nonk_valve
      USE Zmcp         , ONLY: num_mcploc,num_mcpface,n_face_mcp,mcp_on,mapping_mcp,nonk_mcp
      
      IMPLICIT NONE
!
!.....Input
      REAL(8),DIMENSION(ncell_fp) :: s
!.....Output
      REAL(8),DIMENSION(nf_totk,ndim) :: fie_nf
!.....Local variables
      INTEGER :: i,i1,ii,kk,nf_number,istart,zz,tt,i0,istart0,k,nv,ik
      REAL(8) :: fie
!      
!......valve model
!  
     nv=0
     nf_number=nf_number_id(nv)
     istart=istart_nf(1,nf_number)      
     istart0=istart_nfs(0)          
      DO zz=1,num_valveloc                
         tt=mapping_valve(zz)             
         IF(valve_closed(tt).eq.0) CYCLE   
         DO i=1,num_valveface(zz)
            i1=n_face_valve(zz,i) 
            i0=istart0-istart+i1    !i0=istart0+i=istart0+(i1-istart)  
            ii=left_nf(i1)
            fie=s(ii)
            fie_nf(i0,:)=fie*sv_nf(i1,:)             

!           Asymmetric face    
            kk=right_non(i1-istart) !i1=istart+i
            IF(kk.le.ncell_fluid) THEN            
               ik=nonk_valve(zz,i)      ! k=right_nb_k(ik) --> left_nf(k), right_non(k)        
               k=right_nb_k(ik)
               kk=right_non(k)  !right_non(i1-istart)=right_non(k)
               fie=s(kk) 
               fie_nf(ik,:)=-fie*sv_nf(k,:)      
            ENDIF   
         ENDDO
      ENDDO        
!
!......choke flow
!
      IF(choke) THEN
         nv=0
         nf_number=nf_number_id(nv)
         istart=istart_nf(1,nf_number)      
         istart0=istart_nfs(0)          
         DO i=1,num_throatface
            i1=n_face_throat(i)
            i0=istart0-istart+i1    !i0=istart0+i=istart0+(i1-istart) 
            ii=left_nf(i1)
            fie=s(ii) 
            fie_nf(i0,:)=fie*sv_nf(i1,:)
!
!           Asymmetric face
            kk=right_non(i1-istart) !i1=istart+i
            IF(kk.le.ncell_fluid) THEN            
               ik=nonk_throat(i)      ! k=right_nb_k(ik) --> left_nf(k), right_non(k)        
               k=right_nb_k(ik)
               kk=right_non(k)  !right_non(i1-istart)=right_non(k)
               fie=s(kk) 
               fie_nf(ik,:)=-fie*sv_nf(k,:)
            ENDIF   
!
         ENDDO
      ENDIF   
!         
!......MCP model
!     
      nv=0
      nf_number=nf_number_id(nv)
      istart=istart_nf(1,nf_number)      
      istart0=istart_nfs(0)
      DO zz=1,num_mcploc                  !local loop
         tt=mapping_mcp(zz)               !mapping local MCP to global MCP number
         IF(mcp_on(tt).eq.0) CYCLE         !mcp_on: global MCP number    
         DO i=1,num_mcpface(zz)
            i1=n_face_mcp(zz,i) 
            i0=istart0-istart+i1    !i0=istart0+i=istart0+(i1-istart) 
            ii=left_nf(i1)
            fie=s(ii)  
            fie_nf(i0,:)=fie*sv_nf(i1,:)             
!
!           Asymmetric face 
            kk=right_non(i1-istart) !i1=istart+i
            IF(kk.le.ncell_fluid) THEN
               ik=nonk_mcp(zz,i)      ! k=right_nb_k(ik) --> left_nf(k), right_non(k)        
               k=right_nb_k(ik)
               kk=right_non(k)  !right_non(i1-istart)=right_non(k)
               fie=s(kk)
               fie_nf(ik,:)=-fie*sv_nf(k,:)    !k is equal to i1        
            ENDIF   
         ENDDO
      ENDDO    
!         
!
    END SUBROUTINE fluxBC_gradp
    
    SUBROUTINE fluxBC_gradpK1(s,fie_nf)

      USE VOL_DATA     , ONLY: cell   
      USE Zmpi         , ONLY: ncell_fp
      USE Zzone        , ONLY: ncell_fluid
      USE Zparam       , ONLY: ndim    
      USE Zvec_param   , ONLY: nf_totk
      USE Znum_cell    , ONLY: istart_nf,nf_number_id,istart_nfs,right_nb_k
      USE Zvec_index   , ONLY: left_nf,right_non      
      USE Zvec_geo     , ONLY: sv_nf,dxfc_nf,dxfc_non_k
      USE Zrv_choke    , ONLY: num_throatface,n_face_throat,choke,nonk_throat
      USE Zvalve       , ONLY: valve_closed,num_valveloc,num_valveface,n_face_valve,mapping_valve,nonk_valve
      USE Zmcp         , ONLY: num_mcploc,num_mcpface,n_face_mcp,mcp_on,mapping_mcp,nonk_mcp
      USE Zconst2      , ONLY: grav,gfactor
      
      IMPLICIT NONE
!
!.....Input
      REAL(8),DIMENSION(ncell_fp) :: s
!.....Output
      REAL(8),DIMENSION(nf_totk,ndim) :: fie_nf
!.....Local variables
      INTEGER :: i,i1,ii,kk,nf_number,istart,zz,tt,nv,istart0,i0,ik,k
      REAL(8) :: fie,dp
!      
!......valve model
!  
      nv=0
      nf_number=nf_number_id(nv)
      istart=istart_nf(1,nf_number)      
      istart0=istart_nfs(0)        
      DO zz=1,num_valveloc                
         tt=mapping_valve(zz)             
         IF(valve_closed(tt).eq.0) CYCLE   
         DO i=1,num_valveface(zz)
            i1=n_face_valve(zz,i)
            i0=istart0-istart+i1  !i0=istart0+i=istart0+(i1-start) 
            ii=left_nf(i1)
!            
            dp= grav(1)*dxfc_nf(i1,1)*gfactor(ii) &
                +grav(2)*dxfc_nf(i1,2)*gfactor(ii)
            IF(ndim.eq.3) dp=dp+grav(3)*dxfc_nf(i1,3)*gfactor(ii)
            dp=cell%rhom(ii)*dp            
            
            fie=s(ii)+dp
            fie_nf(i0,:)=fie*sv_nf(i1,:)
      
!           Asymmetric face    
            kk=right_non(i1-istart) !i1=istart+i
            IF(kk.le.ncell_fluid) THEN
               ik=nonk_valve(zz,i)      ! k=right_nb_k(ik) --> left_nf(k), right_non(k)        
               k=right_nb_k(ik)
               kk=right_non(k)  !right_non(i1-istart)=right_non(k)
               
               dp= grav(1)*dxfc_non_k(k,1)*gfactor(kk) &
                   +grav(2)*dxfc_non_k(k,2)*gfactor(kk)
               IF(ndim.eq.3) dp=dp+grav(3)*dxfc_non_k(k,3)*gfactor(kk)
               dp=cell%rhom(kk)*dp
            
               fie=s(kk)+dp 
               fie_nf(ik,:)=-fie*sv_nf(k,:)
            ENDIF
      
         ENDDO
      ENDDO        
!
!......choke flow
!
      IF(choke) THEN
         nv=0
         nf_number=nf_number_id(nv)
         istart=istart_nf(1,nf_number)      
         istart0=istart_nfs(0)             
         DO i=1,num_throatface
            i1=n_face_throat(i)
            i0=istart0-istart+i1  !i0=istart0+i=istart0+(i1-start) 
            ii=left_nf(i1)
!            
            dp= grav(1)*dxfc_nf(i1,1)*gfactor(ii) &
                +grav(2)*dxfc_nf(i1,2)*gfactor(ii)
            IF(ndim.eq.3) dp=dp+grav(3)*dxfc_nf(i1,3)*gfactor(ii)
            dp=cell%rhom(ii)*dp            
!            
            fie=s(ii)+dp 
            fie_nf(i0,:)=fie*sv_nf(i1,:)
!
!           Asymmetric face
            kk=right_non(i1-istart) !i1=istart+i
            IF(kk.le.ncell_fluid) THEN
               ik=nonk_throat(i)      ! k=right_nb_k(ik) --> left_nf(k), right_non(k)        
               k=right_nb_k(ik)
               kk=right_non(k)  !right_non(i1-istart)=right_non(k)   
!               
               dp= grav(1)*dxfc_non_k(k,1)*gfactor(kk) &
                   +grav(2)*dxfc_non_k(k,2)*gfactor(kk)
               IF(ndim.eq.3) dp=dp+grav(3)*dxfc_non_k(k,3)*gfactor(kk)
               dp=cell%rhom(kk)*dp            
!            
               fie=s(kk)+dp 
               fie_nf(ik,:)=-fie*sv_nf(k,:)
            ENDIF   
!
         ENDDO
      ENDIF   
!         
!......MCP model
!        
      nv=0
      nf_number=nf_number_id(nv)
      istart=istart_nf(1,nf_number)      
      istart0=istart_nfs(0)      
      DO zz=1,num_mcploc                  !local loop
         tt=mapping_mcp(zz)               !mapping local MCP to global MCP number
         IF(mcp_on(tt).eq.0) CYCLE         !mcp_on: global MCP number    
         DO i=1,num_mcpface(zz)
            i1=n_face_mcp(zz,i) 
            i0=istart0-istart+i1  !i0=istart0+i=istart0+(i1-start) 
            ii=left_nf(i1)
!            
            dp= grav(1)*dxfc_nf(i1,1)*gfactor(ii) &
                +grav(2)*dxfc_nf(i1,2)*gfactor(ii)
            IF(ndim.eq.3) dp=dp+grav(3)*dxfc_nf(i1,3)*gfactor(ii)
            dp=cell%rhom(ii)*dp            
!            
            fie=s(ii)+dp  
            fie_nf(i0,:)=fie*sv_nf(i1,:)             
!
!           Asymmetric face 
            kk=right_non(i1-istart) !i1=istart+i
            IF(kk.le.ncell_fluid) THEN
               ik=nonk_mcp(zz,i)      ! k=right_nb_k(ik) --> left_nf(k), right_non(k)        
               k=right_nb_k(ik)
               kk=right_non(k)  !right_non(i1-istart)=right_non(k)

               dp= grav(1)*dxfc_non_k(k,1)*gfactor(kk) &
                   +grav(2)*dxfc_non_k(k,2)*gfactor(kk)
               IF(ndim.eq.3) dp=dp+grav(3)*dxfc_non_k(k,3)*gfactor(kk)
               dp=cell%rhom(kk)*dp
               
               fie=s(kk)+dp 
               fie_nf(ik,:)=-fie*sv_nf(k,:)    !k is equal to i1        
            ENDIF   
!
         ENDDO
      ENDDO    
!
    END SUBROUTINE fluxBC_gradpK1
!
    SUBROUTINE fluxBC_gradpK2(s,dgdx,fie_nf)

      USE Zmpi         , ONLY: ncell_fp
      USE Zzone        , ONLY: ncell_fluid
      USE Zparam       , ONLY: ndim    
      USE Znum_cell    , ONLY: istart_nf,right_nb_k, &
                               nf_number_id,istart_nfs
      USE Zvec_param   , ONLY: nf_totk
      USE Zvec_geo     , ONLY: sv_nf,dxfc_nf,dxfc_non_k
      USE Zvec_index   , ONLY: left_nf,right_non      
      USE Zrv_choke    , ONLY: num_throatface,n_face_throat,choke,nonk_throat  
      USE Zvalve       , ONLY: valve_closed,num_valveloc,num_valveface,n_face_valve,mapping_valve,nonk_valve
      USE Zmcp         , ONLY: num_mcploc,num_mcpface,n_face_mcp,mcp_on,mapping_mcp,nonk_mcp
      
      IMPLICIT NONE
!
!.....Input
      REAL(8),DIMENSION(ncell_fp) :: s
      REAL(8),DIMENSION(ncell_fp,ndim) :: dgdx
!.....Output
      REAL(8),DIMENSION(nf_totk,ndim) :: fie_nf
!.....Local variables
      INTEGER :: i,i1,ii,kk,nf_number,istart,zz,tt,nv,i0,ik,istart0,k
      REAL(8) :: fie,dp
!      
!......valve model
!  
      nv=0
      nf_number=nf_number_id(nv)
      istart=istart_nf(1,nf_number)      
      istart0=istart_nfs(0)              
      DO zz=1,num_valveloc                
         tt=mapping_valve(zz)             
         IF(valve_closed(tt).eq.0) CYCLE   
         DO i=1,num_valveface(zz)
            i1=n_face_valve(zz,i) 
            i0=istart0-istart+i1  !i0=istart0+i=istart0+(i1-start) 
            ii=left_nf(i1)
!            
            dp= dgdx(ii,1)*dxfc_nf(i1,1) &
               +dgdx(ii,2)*dxfc_nf(i1,2)
            IF(ndim.eq.3) dp=dp+dgdx(ii,3)*dxfc_nf(i1,3)
            fie=s(ii)+dp
            fie_nf(i0,:)=fie*sv_nf(i1,:)            
 
!           Asymmetric face    
            kk=right_non(i1-istart) !i1=istart+i
            IF(kk.le.ncell_fluid) THEN
               ik=nonk_valve(zz,i)      ! k=right_nb_k(ik) --> left_nf(k), right_non(k)        
               k=right_nb_k(ik)
               kk=right_non(k)  !right_non(i1-istart)=right_non(k)
!               
               dp= dgdx(kk,1)*dxfc_non_k(k,1) &
                  +dgdx(kk,2)*dxfc_non_k(k,2)
               IF(ndim.eq.3) dp=dp+dgdx(kk,3)*dxfc_non_k(k,3)
               fie=s(kk)+dp 
               fie_nf(ik,:)=-fie*sv_nf(k,:)
            ENDIF

         ENDDO
      ENDDO        
!
!......choke flow
!
      IF(choke) THEN
         nv=0
         nf_number=nf_number_id(nv)
         istart=istart_nf(1,nf_number)      
         istart0=istart_nfs(0)             
         DO i=1,num_throatface
            i1=n_face_throat(i)
            i0=istart0-istart+i1  !i0=istart0+i=istart0+(i1-start) 
            ii=left_nf(i1)
!            
            dp= dgdx(ii,1)*dxfc_nf(i1,1) &
               +dgdx(ii,2)*dxfc_nf(i1,2)
            IF(ndim.eq.3) dp=dp+dgdx(ii,3)*dxfc_nf(i1,3)          
            fie=s(ii)+dp 
            fie_nf(i0,:)=fie*sv_nf(i1,:)
!
!           Asymmetric face 
            kk=right_non(i1-istart) !i1=istart+i
            IF(kk.le.ncell_fluid) THEN
               ik=nonk_throat(i)      ! k=right_nb_k(ik) --> left_nf(k), right_non(k)        
               k=right_nb_k(ik)
               kk=right_non(k)  !right_non(i1-istart)=right_non(k)
!
               dp= dgdx(kk,1)*dxfc_non_k(k,1) &
                  +dgdx(kk,2)*dxfc_non_k(k,2)
               IF(ndim.eq.3) dp=dp+dgdx(kk,3)*dxfc_non_k(k,3)             
!
               fie=s(kk)+dp 
               fie_nf(ik,:)=-fie*sv_nf(k,:)    !k is equal to i1        
            ENDIF              
!
         ENDDO
      ENDIF   
!         
!......MCP model
!        
      nv=0
      nf_number=nf_number_id(nv)
      istart=istart_nf(1,nf_number)      
      istart0=istart_nfs(0)        
      DO zz=1,num_mcploc                  !local loop
         tt=mapping_mcp(zz)               !mapping local MCP to global MCP number
         IF(mcp_on(tt).eq.0) CYCLE         !mcp_on: global MCP number    
         DO i=1,num_mcpface(zz)
            i1=n_face_mcp(zz,i) 
            i0=istart0-istart+i1  !i0=istart0+i=istart0+(i1-start) 
            ii=left_nf(i1)
!            
            dp= dgdx(ii,1)*dxfc_nf(i1,1) &
               +dgdx(ii,2)*dxfc_nf(i1,2)
            IF(ndim.eq.3) dp=dp+dgdx(ii,3)*dxfc_nf(i1,3)  
            fie=s(ii)+dp  
            fie_nf(i0,:)=fie*sv_nf(i1,:)             

!           Asymmetric face 
            kk=right_non(i1-istart) !i1=istart+i
            IF(kk.le.ncell_fluid) THEN
               ik=nonk_mcp(zz,i)      ! k=right_nb_k(ik) --> left_nf(k), right_non(k)        
               k=right_nb_k(ik)
               kk=right_non(k)  !right_non(i1-istart)=right_non(k)
!
               dp= dgdx(kk,1)*dxfc_non_k(k,1) &
                  +dgdx(kk,2)*dxfc_non_k(k,2)
               IF(ndim.eq.3) dp=dp+dgdx(kk,3)*dxfc_non_k(k,3)             
!
               fie=s(kk)+dp 
               fie_nf(ik,:)=-fie*sv_nf(k,:)    !k is equal to i1        
            ENDIF   
!
         ENDDO
      ENDDO
!
      END SUBROUTINE fluxBC_gradpK2    
