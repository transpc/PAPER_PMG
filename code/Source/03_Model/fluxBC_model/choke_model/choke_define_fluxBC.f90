!
      SUBROUTINE choke_define_fluxBC
!
!     DATA STRUCTURE
!
!     num_throatface=total number of faces at a broken throat
!     fzone_throat(1,1)=fluid zone index inside the throat (rv_parameters.in)
!     fzone_throat(1,2)=fluid zone index outside the throat (rv_parameters.in)
!     icell_throat(1:num_thratface)=cell index just inside the throat
!     ocell_throat(1:num_thratface)=cell index just outside the throat
!     n_face_throat(1:num_thratface)=non-face index of the throat     
!     dir_face_throat(1:num_thratface)=non face direction of the throat
!     ick_dir=choking flow direction
!     theta_avg=choking face direction (rad.) (90deg=facing upward, -90deg=facing downward, 0deg=horizontal)      
!  
      USE Znum_cell    , ONLY: istart_nf  
      USE Zvec_index   , ONLY: left_nf,right_non
      USE Zvec_geo     , ONLY: saa_nf
      USE Zvec_param   , ONLY: nf_nonk
      USE Znum_cell    , ONLY: right_nb_k      
      USE Zrv_choke    
      USE Zzone        , ONLY: nzone
      USE Zcore        , ONLY: np,myrank
      USE Zcoord1      , ONLY: xloc
      USE Zparam       , ONLY: pi
      USE Zinterface      
      USE Zmpi
      USE Zzone        , only: ncell_fluid_all
!
      IMPLICIT NONE
!
      INTEGER :: get_global_cell             
      INTEGER :: nf_number,istart,isize,i1,ir,j,k,nn
      INTEGER :: i,ii,kk,ii_throat,kk_throat,ii_tmp,kk_tmp
      INTEGER :: ithroat,othroat    
      INTEGER :: err,unit_valve,num_zone,iig,kkg,k1,k2
      INTEGER,ALLOCATABLE ::cell_fluxBC(:),zone_fluxBC(:),nzone_fluxBC(:)      
      LOGICAL,SAVE :: initial=.true.
      LOGICAL,SAVE :: initial_nzone=.true.      
!
!.....Open FluxBCzone file
!      
      unit_valve=423 
      OPEN(unit_valve,file='FluxBCZone.in',status='old',iostat=err)
      IF(err.eq.0)then
         IF(myrank.eq.0)WRITE(*       ,"(11x,a)")'Reading FluxBCZone_valve.in...'
         READ(unit_valve,*) num_zone
         ALLOCATE(cell_fluxBC(num_zone),zone_fluxBC(num_zone),nzone_fluxBC(ncell_fluid_all))
         nzone_fluxBC=0
         DO k=1,num_zone
            READ(unit_valve,*) cell_fluxBC(k),zone_fluxBC(k)
            ii=cell_fluxBC(k)
            nzone_fluxBC(ii)=zone_fluxBC(k) !ii=global, nzone_valve=FluxBC zone number
         ENDDO
      ELSE
         print*,'>> ERROR! FluxBCZone.in file must be in the folder.'
         pause
         stop
      ENDIF
      CLOSE(unit_valve)      
!
!......Define global nzone
!      
      IF(initial_nzone) THEN
         PRINT*,'>>>Choking face is finding...'
         IF(np.gt.1) CALL communicate_1d_int(nzone)         
         initial_nzone=.false.
      ENDIF         
!      
!.....Initialize the total number of throat face
!
      num_throatface=0
      nf_number=0
      istart=istart_nf(1,nf_number)
      isize =istart_nf(2,nf_number)
      ithroat=fzone_throat(1,1)     !inside fluid zone number of throat 
      othroat=fzone_throat(1,2)     !outside fluid zone number of throat 
      DO i=1,isize
         i1=istart+i
         ii=left_nf(i1)
         kk=right_non(i)
         iig=get_global_cell(ii) !global
         kkg=get_global_cell(kk) !global            
!               
         k1=0
         k2=0
         DO k=1,num_zone
            IF(iig.eq.cell_fluxBC(k)) THEN
               k1=k
            ENDIF
            IF(kkg.eq.cell_fluxBC(k)) THEN
               k2=k
            ENDIF               
         ENDDO
         IF(k1.ne.0.and.k2.ne.0) THEN
            IF((zone_fluxBC(k1).eq.ithroat.and.zone_fluxBC(k2).eq.othroat).or. &
               (zone_fluxBC(k2).eq.ithroat.and.zone_fluxBC(k1).eq.othroat)) THEN            
               num_throatface=num_throatface+1               
            ENDIF
         ENDIF          
     ENDDO     
!
!....Initialize the array
!      
      IF(initial) THEN
         ALLOCATE(icell_throat(num_throatface),ocell_throat(num_throatface))
         ALLOCATE(n_face_throat(num_throatface),dir_face_throat(num_throatface))
         ALLOCATE(fluxl_throatface(num_throatface),fluxg_throatface(num_throatface),fluxd_throatface(num_throatface))
         ALLOCATE(nonk_throat(num_throatface))        !local domain          
         initial=.false.
      ENDIF            
!      
!.....Finding throat face between computing cell and pressure outlet cell
!
      num_throatface=0
      nn=0
      nf_number=0
      istart=istart_nf(1,nf_number)
      isize =istart_nf(2,nf_number)          
      ithroat=fzone_throat(1,1)     !inside fluid zone number of throat 
      othroat=fzone_throat(1,2)     !outside fluid zone number of throat          
      DO i=1,isize
         i1=istart+i
         ii=left_nf(i1)   !local
         kk=right_non(i)  !local 
         iig=get_global_cell(ii) !global
         kkg=get_global_cell(kk) !global

         k1=0
         k2=0
         DO k=1,num_zone
            IF(iig.eq.cell_fluxBC(k)) THEN
               k1=k
            ENDIF
            IF(kkg.eq.cell_fluxBC(k)) THEN
               k2=k
            ENDIF               
         ENDDO
         IF(k1.ne.0.and.k2.ne.0) THEN
            IF((zone_fluxBC(k1).eq.ithroat.and.zone_fluxBC(k2).eq.othroat).or. &
               (zone_fluxBC(k2).eq.ithroat.and.zone_fluxBC(k1).eq.othroat)) THEN
!
                num_throatface=num_throatface+1 
                n_face_throat(num_throatface)=i1                 !non-face index for the valve faces                         
! 
               IF(zone_fluxBC(k1).eq.ithroat.and.zone_fluxBC(k2).eq.othroat) THEN
                  dir_face_throat(num_throatface)=1              !owner of the throat face is left
                  icell_throat(num_throatface)=ii                !fzone_throat(1) local ii
                  ocell_throat(num_throatface)=kk                !fzone_throat(2) local kk (kk=1~ncell_fp)

                  IF(xloc(ii,1).ne.xloc(kk,1)) THEN
                      ick_dir=1
                      theta_avg=0.d0
                  ELSEIF(xloc(ii,2).ne.xloc(kk,2)) THEN 
                      ick_dir=2
                      theta_avg=0.d0
                  ELSEIF(xloc(ii,3).le.xloc(kk,3)) THEN 
                      ick_dir=3
                      theta_avg=pi/2.d0
                  ELSEIF(xloc(ii,3).ge.xloc(kk,3)) THEN 
                      ick_dir=3
                      theta_avg=-pi/2.d0                   
                  ENDIF                  
                  
               ELSEIF(zone_fluxBC(k2).eq.ithroat.and.zone_fluxBC(k1).eq.othroat) THEN
                  dir_face_throat(num_throatface)=-1
                  icell_throat(num_throatface)=kk               !fzone_throat(1) local kk (kk=1~ncell_fp)
                  ocell_throat(num_throatface)=ii               !fzone_throat(2) local ii

                  IF(xloc(kk,1).ne.xloc(ii,1)) THEN
                      ick_dir=1
                      theta_avg=0.d0
                  ELSEIF(xloc(kk,2).ne.xloc(ii,2)) THEN 
                      ick_dir=2
                      theta_avg=0.d0
                  ELSEIF(xloc(kk,3).le.xloc(ii,3)) THEN 
                      ick_dir=3
                      theta_avg=pi/2.d0
                  ELSEIF(xloc(kk,3).ge.xloc(ii,3)) THEN 
                      ick_dir=3
                      theta_avg=-pi/2.d0                   
                  ENDIF                  
                  
               ENDIF                      
!                     
            ENDIF
         ENDIF
      ENDDO
!      
      IF(num_throatface.gt.1000) THEN
         PRINT*,'ERROR!'
         PRINT*,'throatface number is over the pre-defined array size (maxface).'
         PAUSE
         STOP
      ENDIF          
!      
!.....Calculate the throat area
!      
      throat_area=0.d0
      DO i=1,num_throatface
          i1=n_face_throat(i)
          throat_area=throat_area+saa_nf(i1)
      ENDDO
!      
!.....Find nonk_throat for asymmetric sumation
!            
      nf_number=0
      istart=istart_nf(1,nf_number)  
      DO i=1,num_throatface
         i1=n_face_throat(i) 
         ir=i1-istart    !i1=istart+i
         ii=left_nf(i1)
         kk=right_non(ir)
!            
         ii_throat=icell_throat(i)
         kk_throat=ocell_throat(i)
!            
         DO j=1,nf_nonk
            k=right_nb_k(j) 
            ii_tmp=left_nf(k)               
            kk_tmp=right_non(k) 
!               
            IF((ii_throat.eq.ii_tmp.and.kk_throat.eq.kk_tmp).or.(ii_throat.eq.kk_tmp.and.kk_throat.eq.ii_tmp)) THEN
               nonk_throat(i)=j  
            ENDIF
!                   
         ENDDO              
      ENDDO
!      
      DEALLOCATE(cell_fluxBC,zone_fluxBC,nzone_fluxBC)      
!      
      RETURN
      END SUBROUTINE choke_define_fluxBC
    
!
      SUBROUTINE choke_define_fluxBC_apr1400
!
!     DATA STRUCTURE
!
!     num_throatface=total number of faces at a broken throat
!     fzone_throat(1,1)=fluid zone index inside the throat (rv_parameters.in)
!     fzone_throat(1,2)=fluid zone index outside the throat (rv_parameters.in)
!     icell_throat(1:num_thratface)=cell index just inside the throat
!     ocell_throat(1:num_thratface)=cell index just outside the throat
!     n_face_throat(1:num_thratface)=non-face index of the throat     
!     dir_face_throat(1:num_thratface)=non face direction of the throat
!     ick_dir=choking flow direction
!     theta_avg=choking face direction (rad.) (90deg=facing upward, -90deg=facing downward, 0deg=horizontal)      
!  
      USE Znum_cell    , ONLY: istart_nf  
      USE Zvec_index   , ONLY: left_nf,right_non
      USE Zvec_geo     , ONLY: saa_nf
      USE Zbc_index    , ONLY: npb 
      USE Zrv_choke    
      USE Zparam       , ONLY: pi
      USE Zinterface      
      USE Zmpi
!
      IMPLICIT NONE
!
      INTEGER :: nf_number,istart,isize,i1
      INTEGER :: i,ii,kk
      LOGICAL,SAVE :: initial=.true.
!      
!.....Initialize the total number of throat face
!
      num_throatface=0
      nf_number=0
      istart=istart_nf(1,nf_number)
      isize =istart_nf(2,nf_number)
      DO i=1,isize
         i1=istart+i
         ii=left_nf(i1)
         kk=right_non(i)
         IF((npb(kk).eq.1.and.npb(ii).ne.1).or.(npb(kk).ne.1.and.npb(ii).eq.1)) THEN
            num_throatface=num_throatface+1
         ENDIF
      ENDDO          
!
!....Initialize the array
!      
      IF(initial) THEN
         ALLOCATE(icell_throat(num_throatface))
         ALLOCATE(n_face_throat(num_throatface),dir_face_throat(num_throatface))
         ALLOCATE(fluxl_throatface(num_throatface),fluxg_throatface(num_throatface),fluxd_throatface(num_throatface))
         initial=.false.
      ENDIF            
!      
      IF(num_throatface.gt.1000) THEN
         PRINT*,'ERROR!'
         PRINT*,'throatface number is over the pre-defined array size (maxface).'
         PAUSE
         STOP
      ENDIF          
!      
!.....Special case using pboun for choking condition
!      
      num_throatface=0
      nf_number=0
      istart=istart_nf(1,nf_number)
      isize =istart_nf(2,nf_number)
      DO i=1,isize
         i1=istart+i
         ii=left_nf(i1)
         kk=right_non(i)
         IF((npb(kk).eq.1.and.npb(ii).ne.1).or.(npb(kk).ne.1.and.npb(ii).eq.1)) THEN
            num_throatface=num_throatface+1
            n_face_throat(num_throatface)=i1
            IF((npb(kk).eq.1.and.npb(ii).ne.1))dir_face_throat(num_throatface)=1
            IF((npb(kk).ne.1.and.npb(ii).eq.1))dir_face_throat(num_throatface)=-1
         ENDIF
         IF(npb(kk).eq.1.and.npb(ii).ne.1) THEN
            icell_throat(num_throatface)=ii                    
         ENDIF  
      ENDDO          
      ALLOCATE(ul_throatface(num_throatface),ug_throatface(num_throatface),ud_throatface(num_throatface))
!      
!.....Calculate the throat area
!      
      throat_area=0.d0
      DO i=1,num_throatface
          i1=n_face_throat(i)
          throat_area=throat_area+saa_nf(i1)
      ENDDO      
!      
      RETURN
    END SUBROUTINE choke_define_fluxBC_apr1400    
    
    
    
    
 