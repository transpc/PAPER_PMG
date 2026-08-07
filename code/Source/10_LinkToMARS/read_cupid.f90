!
     SUBROUTINE read_cupid(time_system)
!                                                                       
!     Processes the CUPID input data
!                                                                                                                   
      USE Zcore        , only:myrank,np,cupid_mars,cupid_alone
      USE Zconst1      , ONLY:cplmars,cplmaster
      USE Zconst2      , only:i_repeat
      USE Ztimecon     , only:time      
      USE c3com_cupid  , only:i3cupid
      USE Zvec_param   , ONLY: vl_f
      USE Zzone        , ONLY: ncell_cond_all
      USE Zmars        , ONLY: time_mars
      USE Zio_unit     , ONLY: unit_somaflow,unit_log
!
      IMPLICIT NONE
!DEC$IF defined (MCC)      
      INCLUDE 'c3com.h'
!DEC$ELSEIF defined (MCC_DLL)      
      INCLUDE 'c3com.h'
      !dec$ attributes dllexport :: read_cupid      
!DEC$ELSEIF defined (SPACE)          
      INCLUDE 'c3com.h'
      INCLUDE 'c3com_space.h'
      !dec$ attributes dllexport :: read_cupid
!DEC$ENDIF          
!     
      INTEGER i,nout,tripin   
      REAL(8) time_system
!
!.....Get important parameter
!
      myrank=i3myrank
      np=i3np
      cupid_alone=i3cupid_alone
      cupid_mars=i3cupid_mars
!
!...Check whether the 3D data exists
!
!DEC$IF defined (MCC)     
      IF(myrank.eq.0)open(662,file='MARS_volume_index.dat')
      IF(i3nic(2).eq.0)THEN
         cplmars=0
      i3cplmars=cplmars
!         CALL read_cupid_mars
         RETURN
      ENDIF   
!DEC$ELSEIF defined (MCC_DLL)      
      IF(myrank.eq.0)open(662,file='MARS_volume_index.dat')
      IF(i3nic(2).eq.0)THEN
         cplmars=0
      i3cplmars=cplmars
!         CALL read_cupid_mars
         RETURN
      ENDIF        
!DEC$ELSEIF defined (SPACE) 
      cupid_mars=1 !cupid_space_debug
      np=1         !cupid_space_debug
      myrank=0     !cupid_space_debug 
      CALL c3com_copy_S2C
!DEC$ENDIF          
      CALL barrier_mpi
!
!...open input and output files
!
      i_repeat=0
      CALL open_files 
!
!.....input data and initialization
!
      nout=-1
      CALL read_flow(nout)  
!      
      CALL read_grid
!
!to enable vector code
      CALL vectorize_index
      IF(ncell_cond_all.gt.0) CALL vectorize_index_solid
      CALL vectorize_index_connectivity
      CALL vectorize_geo_var
      CALL vectorize_major_flux
!
      nout=0 
      CALL read_flow(nout)!  
      CLOSE(unit_somaflow) 
!
      CALL vectorize_scalar_upwind
!      
      CALL vectorize_deallocate_face      
!
!DEC$IF defined (MCC)      
      time=c3time_sys
      time_mars=c3time_sys  
      IF(myrank.eq.0)WRITE(*,"(10x,a,1e12.5,a,1e12.5)")'cupid_time=',time,'mars_time=',time_mars
!DEC$ELSEIF defined (MCC_DLL)          
      time=c3time_sys
      time_mars=c3time_sys  
      IF(myrank.eq.0)WRITE(*,"(10x,a,1e12.5,a,1e12.5)")'cupid_time=',time,'mars_time=',time_mars      
!DEC$ELSEIF defined (SPACE)          
      CALL broadcast_r(time_system,1) !scc-mpi
       time=time_system
      IF(myrank.eq.0)WRITE(*,"(10x,a,1e12.5,a,1e12.5)")'cupid_time=',time,'SPACE_time=',time_system
      IF(i3nic(2).gt.0)CALL set_cupid_cell_no   !MCC-jjj
!DEC$ENDIF 
!
!      CALL read_cupid_mars
!
      i3cplmaster=cplmaster  
      i3cplmars=cplmars
!
      RETURN 
      END SUBROUTINE read_cupid
!------------------------------------------------------------------------------------------------
!DEC$IF defined (SPACE)      
!jhhong - space-cupid coupling - for sharing of common block c3com
     subroutine return_addr(addr)
     INCLUDE 'c3com_space.h'
!dec$ attributes dllexport :: return_addr
     integer addr
!   cupid --> space 로 옮겨줘야 할 값
!    
     addr = LOC(s3dt_super)
!     
     end subroutine return_addr
!DEC$ENDIF
!------------------------------------------------------------------------------------------------
