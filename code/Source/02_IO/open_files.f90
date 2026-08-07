!
      SUBROUTINE open_files
!
!     This routine opens input and ouptut files
!
      USE Zcore       , ONLY: cupid_alone,cupid_mars,np,myrank
      USE Zparam      , ONLY: ndim,nn,ns,outfilename,mesh_openfoam,mesh_binary
      USE Zconst1     , ONLY: restart,vv_prob
      USE Zio_unit    , ONLY: unit_somaflow,unit_grid,unit_saveout,unit_log,unit_screen
      !,mboron,mdiffoff,cplmars,cplmaster,topolsurface,mdiffscheme
      USE Zmars       , ONLY: ncupvol,ncell_old
      USE Ztimecon    , ONLY: itim_last,nbline !,iso_thermal
      USE Zconst1     , ONLY: restart
      USE viewData_common , ONLY: viewField      
      !USE Zgradoption , ONLY: ifrink
      USE Znode        , ONLY: nmax_vertex,n_node
!      USE Zporous     , ONLY: vfporous
      USE unitManager  , ONLY: createUnit
!
      IMPLICIT NONE
!      
      CHARACTER(8) :: date
      CHARACTER(10):: times
      CHARACTER(5) :: zone
!      
      INTEGER,DIMENSION(8) :: values 
      INTEGER i,IX,ios
      INTEGER npgrid,nmesh
      INTEGER flowin,gridin,savein     
      INTEGER foam1,n_face,n_bface,n_bc,n_zone
      INTEGER num_fzone_count,num_szone_count
      INTEGER :: itmp(9)
!            
      INTEGER save_option,noutput,iprn,itim
      REAL(8) tplot_dt, tplot_num, tplot_prop
      NAMELIST /problem_description/ vv_prob, ndim, num_fzone_count, num_szone_count
      NAMELIST /post/ save_option,noutput,iprn,restart,itim,tplot_dt,tplot_num,tplot_prop,viewField
!      NAMELIST /misc_option/ mboron,mdiffoff,cplmars,cplmaster,topolsurface,mdiffscheme,iso_thermal,ifrink,vfporous,&
!                              vfporous,i_droplet,i_fs_temp_intpol
         
!      
      REAL dum(4)
!
      DATA flowin,gridin/1,1/
      itim_last=0
!      
!.....List of files
!     (Input)
!     4  : somaGrid.in
!     5  : somaFlow.in 
!     16 : tpfh2o
!
!     (history (time, dt, mass error))
!     97 : log.dat
!     
!     (Saveout indicator)
!     222: saveout.dat
!     44 :restart.dat 
!
!     (User-defined) 
!     101:inlet_bcdat.inp, 60:rocom_transient.dat
!     59:cboron.dat, 51~56:height01~06, 61~66:Nheight01~06, 67~68:Nradius01~02
!     20:dam_front.dat or sbloca.dat, 30:plot.out, 50:tec2d.dat, 111:sa.dat, 333:ta.dat or osc.dat
!
!     (MASTER coupling)
!     2010:3dkin.n
!
!.....Open 'somaFlow.in'
!
      unit_somaflow=createUnit("somaFlow.in")
      unit_somaflow=812
      OPEN(unit_somaflow, file='somaFlow.in', recl=200, delim='APOSTROPHE', iostat=flowin)
      IF(flowin.ne.0)then
         WRITE(*,*)'Program was terminated due to lack of <somaFlow.in>!'
         WRITE(unit_log,*)'Program was terminated due to lack of <somaFlow.in>!'
         STOP
      ENDIF
      READ(unit_somaflow, nml=problem_description)
!      READ(812, nml=misc_option)      
      READ(unit_somaflow, nml=post)      
      REWIND(unit_somaflow)
      ix=ndim
!
!.....Open 'somaGrid.in'
!
      mesh_openfoam=0
      mesh_binary=0
!     
      unit_grid=createUnit('grid.in')
      unit_grid=4   
!      
      IF(cupid_alone.eq.1)then
         IF(myrank.eq.0)THEN
            OPEN(unit_grid,file='somaGrid.in',form='binary',status='old',iostat=gridin)
            mesh_binary=1
         ENDIF 
         IF(np.gt.1) CALL broadcast_i1(gridin)
         IF(gridin.ne.0)then
            IF(myrank.eq.0) OPEN(unit_grid,file='foamGrid.in',status='old',iostat=foam1)
            IF(np.gt.1) CALL broadcast_i1(foam1)
            IF(foam1.ne.0)then
               IF(myrank.eq.0) WRITE(*,*)'Program was terminated due to lack of <foamGrid.in>!'
               CALL finalize_mpi
               STOP               
            ELSE
               mesh_openfoam=1
            ENDIF   
            IF(myrank.eq.0) THEN
               WRITE(*,"(11x,a)")'Reading foamGrid.in...'
               READ(unit_grid,*)outfilename
               READ(unit_grid,*)n_node,nn,n_face,n_bface,n_bc,n_zone,nmax_vertex,ns,ndim
               REWIND(unit_grid)
               IF(np.gt.1) THEN
                  itmp(1)=n_node
                  itmp(2)=nn
                  itmp(3)=n_face
                  itmp(4)=n_bface
                  itmp(5)=n_bc
                  itmp(6)=n_zone
                  itmp(7)=nmax_vertex
                  itmp(8)=ns
                  itmp(9)=ndim
               ENDIF
            ENDIF   
            IF(np.gt.1) THEN
               CALL broadcast_i(itmp,9)
               n_node      =itmp(1)
               nn          =itmp(2)
               n_face      =itmp(3)
               n_bface     =itmp(4)
               n_bc        =itmp(5)
               n_zone      =itmp(6)
               nmax_vertex =itmp(7)
               ns          =itmp(8)
               ndim        =itmp(9)
            ENDIF   
         ELSE
            IF(myrank.eq.0) THEN
               WRITE(*,"(11x,a)")'Reading somaGrid.in...'
               READ(unit_grid)outfilename
               READ(unit_grid)nn,ndim,ns
               itmp(1)=nn
               itmp(2)=ndim
               itmp(3)=ns
            ENDIF   
            IF(np.gt.1) THEN
               CALL broadcast_i(itmp,3)
               nn  =itmp(1)
               ndim=itmp(2)
               ns  =itmp(3)
            ENDIF   
            ncell_old(:)=nn            
         ENDIF         
      ELSEIF(cupid_mars.eq.1)then
         OPEN(unit_grid,file='somaGrid9to1mc.in',form='binary',status='old',iostat=gridin)
         mesh_binary=1
         IF(gridin.ne.0)then
            mesh_binary=0
            OPEN(unit_grid,file='somaGrid.in',status='old',iostat=gridin)
            IF(gridin.ne.0)THEN
               OPEN(unit_grid,file='foamGrid.in',status='old',iostat=foam1)
               IF(foam1.eq.0)then
                  IF(myrank.eq.0)WRITE(*,"(11x,a)")'<foamGrid.in> is used instead of <somaGrid9to1mc.in>!'
                  mesh_openfoam=1
                  READ(unit_grid,*)outfilename
                  READ(unit_grid,*)n_node,nn,n_face,n_bface,n_bc,n_zone,nmax_vertex,ns,ndim
                  REWIND(unit_grid)  
               ELSE             
                  IF(myrank.eq.0)WRITE(*,"(11x,a)")'Program was terminated due to lack of <somaGrid9to1mc.in>!'
                  STOP
               ENDIF          
            ELSE   
               IF(myrank.eq.0)WRITE(*,"(11x,a)")'<SomaGrid.in> is used instead of <somaGrid9to1mc.in>!'
               READ(unit_grid,*)outfilename
               READ(unit_grid,*)nn,ndim,ns !udfn_grid_rearrange               
            ENDIF
         ELSE
            READ(unit_grid)outfilename
            READ(unit_grid)nn,ndim,ns,npgrid,nmesh,ncupvol,ncell_old(1),ncell_old(2),ncell_old(3)  
            IF(ncupvol.le.0)then
               WRITE(*,*)'No cupvols for cupid/mars run!'
               STOP
            ENDIF      
         ENDIF  
      ENDIF  
!
!.....Check dimension
!
      IF(ndim.ne.ix)then
         PRINT *,"Dimension mismatch between T/H and grid inputs(somaFlow.in and somaGriD.in)!"
         PRINT *,"Dimensions in T/H and grid inputs:",ndim,ix
         PAUSE
         STOP
      ENDIF       
      IF(ndim.eq.2)THEN
      ELSEIF(ndim.eq.3)THEN
      ELSE
         PRINT *,"ERROR: number of dimensions should be 2 or 3 !"
         PAUSE
         STOP
      ENDIF            
!
!.....Open output files to write the mass balance, calculation time, time history, event
! 
      restart=mod(restart,10)
      IF(myrank.eq.0) THEN     
         unit_saveout=createUnit("saveout.dat")
         unit_log=createUnit("log.dat")
         unit_saveout=96
         unit_log=97
         IF(restart.eq.0)THEN
            OPEN(unit_saveout,file="saveout.dat")        
            OPEN(unit_log,file="log.dat")             
         ELSE 
            OPEN(unit_saveout,file="saveout.dat",status='old',iostat=savein)
            IF(savein.ne.0)then
               WRITE(*,*)'Program was terminated due to lack of <saveout.dat>!'
               WRITE(unit_log,*)'Program was terminated due to lack of <saveout.dat>!'
               STOP
            ENDIF    
!                   
            DO WHILE(.true.)
               READ(unit_saveout,*,iostat=ios)
               IF(ios.gt.0)THEN
                  WRITE(*,*)'Program was terminated due to the problem in <saveout.dat>!'
                  WRITE(unit_log,*)'Program was terminated due to the problem in <saveout.dat>!'
                  STOP  
               ELSEIF(ios.lt.0)THEN
                  EXIT
               ELSE
                  nbline=nbline+1
               ENDIF
            ENDDO
            REWIND(unit_saveout)
            DO i=1, nbline-1
               READ(unit_saveout,*)
            ENDDO
            READ(unit_saveout,*) dum(1),dum(2),dum(3),itim_last    
            CLOSE(unit_saveout)
!

            OPEN(unit_log,file="log.dat",position='append') 
         ENDIF
!
         OPEN(unit=unit_screen,CARRIAGECONTROL='FORTRAN')
!         
         CALL date_and_time(date,times,zone,values)
         WRITE(unit_log,*) '## Calculation starts at ',date,'(YYYYMMDD), ',times,'(HHMMSS.SSS)'   
!      
         IF(restart.eq.1)THEN
            WRITE(*,*)'          This is a restart run of CUPID.'
           WRITE(unit_log,*)'          This is a restart run of CUPID.'
         ELSEIF(restart.eq.2)THEN
            WRITE(*,*)'          This is a auto restart of CUPID.'
           WRITE(unit_log,*)'          This is a auto restart of CUPID.'
         ELSE
            restart=0
            WRITE(*,*)'          This is a new run of CUPID.'
           WRITE(unit_log,*)'          This is a new run of CUPID.'
         ENDIF
!         
      ENDIF

! 
!.....broadcast itim_last 
!     
      IF(np.gt.1) CALL broadcast_i1(itim_last)
!
      END SUBROUTINE open_files
