      SUBROUTINE GID_out_cell_lsj(time,nout)
!
!.....Save Output for 3D a graphic processor using Open-GL library
!
      USE VOL_DATA    , ONLY: cell
      USE Zmpi        , ONLY: ncell_fp
      USE Zzone       , ONLY: ncell_fluid_all,ncell_fluid
      USE Zparam      , ONLY: ndim
      USE Zcore       , ONLY: myrank
      USE Zconst1     , ONLY: restart
      USE Zpress      , ONLY: p
      USE Zvector     , ONLY: vl_n
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i,j,na,nout
      LOGICAL,SAVE :: initial_restart=.true.
      LOGICAL,SAVE :: initial_start=.true.
      CHARACTER(16) :: meshtype   
      REAL(8) :: time 
!.....Local arrays
      INTEGER,DIMENSION(:),ALLOCATABLE :: regime_global
      REAL(8),DIMENSION(:),ALLOCATABLE :: pr,ag,tl
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: vlp
!
!.....GID OUTPUT PRINT: flavis.res
!
      IF(initial_start.and.ndim.eq.2)THEN
         meshtype='Quadrilateral'
      ENDIF
      IF(initial_start.and.ndim.eq.3)THEN
         meshtype='Hexahedra'
      ENDIF
!
      IF(myrank.eq.0 .and. nout.eq.0)THEN
         WRITE(41,*) 'GiD Post Results File 1.0'
         WRITE(41,*) ''
         WRITE(41,*) 'GaussPoints "Cell Center" ElemType   ',meshtype
         WRITE(41,*) 'Number Of Gauss Points: 1'
         WRITE(41,*) 'Natural Coordinates: Internal'
         WRITE(41,*) 'End gausspoints'
         WRITE(41,*) ''
         WRITE(41,*) ''
         initial_start=.false.
      ENDIF
      IF(myrank.eq.0 .and. restart.ne.0 .and. initial_restart)THEN
         WRITE(41,*) 'GiD Post Results File 1.0'
         WRITE(41,*) ''
         WRITE(41,*) 'GaussPoints "Cell Center" ElemType   ',meshtype
         WRITE(41,*) 'Number Of Gauss Points: 1'
         WRITE(41,*) 'Natural Coordinates: Internal'
         WRITE(41,*) 'End gausspoints'
         WRITE(41,*) ''
         WRITE(41,*) ''
         initial_start=.false.
         initial_restart=.false.
      ENDIF  
!
      na=ncell_fluid_all
      IF(myrank.eq.0) THEN
         ALLOCATE(pr(na),ag(na))
      ELSE
         ALLOCATE(pr(1),ag(1))
      ENDIF
      CALL gatherv_r(p          ,ncell_fluid,pr,na,0)
      CALL gatherv_r(cell%alphag,ncell_fluid,ag,na,0)
!
      IF(myrank.eq.0)then    
!
!........Pressure
!
         WRITE(41,20) time
         WRITE(41,*) 'ComponentNames "P"'
         WRITE(41,*) 'Values'
         DO i=1,ncell_fluid_all
            WRITE(41,12) i,pr(i)
         ENDDO
         WRITE(41,*) 'End Values'
         WRITE(41,*) ''
!         
!........Gas Fraction
!
         WRITE(41,30) time
         WRITE(41,*) 'ComponentNames "Alphag"'
         WRITE(41,*) 'Values'
         DO i=1,ncell_fluid_all
            WRITE(41,12) i,ag(i)
         ENDDO
         WRITE(41,*) 'End Values'
         WRITE(41,*) ''
      ENDIF
      DEALLOCATE(pr,ag) 
!
      IF(myrank.eq.0) THEN
         ALLOCATE(tl(na))
      ELSE
         ALLOCATE(tl(1))
      ENDIF
      CALL gatherv_r(cell%tl,ncell_fluid,tl,ncell_fluid_all,0)
!
      IF(myrank.eq.0)then 
!      
!........Liquid Temperature
!
         WRITE(41,60) time
         WRITE(41,*) 'ComponentNames "Tl"'
         WRITE(41,*) 'Values'
         DO i=1,ncell_fluid_all
            WRITE(41,12) i,tl(i)
         ENDDO
         WRITE(41,*) 'End Values'
         WRITE(41,*) ''
!
      ENDIF
      DEALLOCATE(tl)         
!
!.....Liquid Velocity
!
      IF(myrank.eq.0)then  
         ALLOCATE(vlp(ncell_fluid_all,ndim))
      ELSE
         ALLOCATE(vlp(1,1))
      ENDIF
      CALL gatherv_r_2d(vl_n,ncell_fp,vlp,ncell_fluid,ncell_fluid_all,0)
!
      IF(myrank.eq.0)then
         WRITE(41,80) time
         WRITE(41,*) 'ComponentNames "X-component", "Y-component", "Z-component", "All"'
         WRITE(41,*)'Values'
         DO i=1,ncell_fluid_all
            WRITE(41,12) i,(vlp(i,j),j=1,ndim)
         ENDDO
         WRITE(41,*)'End Values'
         WRITE(41,*)''
      ENDIF   
      DEALLOCATE(vlp)  
!
!.....Regime
!
      IF(myrank.eq.0)then  
         ALLOCATE(regime_global(ncell_fluid_all))
      ELSE
         ALLOCATE(regime_global(1))
      ENDIF
      CALL gatherv_i(cell%regime,ncell_fluid,regime_global,ncell_fluid_all,0)
!
      IF(myrank.eq.0)then  
         WRITE(41,120)time
         WRITE(41,*)'ComponentNames "FlowRegime"'
         WRITE(41,*)'Values'
         DO i=1,ncell_fluid_all
            IF(regime_global(i).eq.11) WRITE(41,*) i,1
            IF(regime_global(i).eq.12) WRITE(41,*) i,2
            IF(regime_global(i).eq.13) WRITE(41,*) i,3
            IF(regime_global(i).eq.21) WRITE(41,*) i,4
            IF(regime_global(i).eq.22) WRITE(41,*) i,5
            IF(regime_global(i).eq.23) WRITE(41,*) i,6
            IF(regime_global(i).eq. 3) WRITE(41,*) i,7
         ENDDO
         WRITE(41,*)'End Values'
         WRITE(41,*)''
      ENDIF     
      DEALLOCATE(regime_global) 
!
   12 FORMAT(i5,5x,7(e20.10,1x))
   20 FORMAT('Result "Pressure" "Pressure"', f20.10, 2x, 'Scalar OnGaussPoints "Cell Center"')
   30 FORMAT('Result "Gas Fraction" "Volume Fraction"', f20.10, 2x, 'Scalar OnGaussPoints "Cell Center"')
   60 FORMAT('Result "Tl" "Temperature"', f20.10, 2x, 'Scalar OnGaussPoints "Cell Center"')
   80 FORMAT('Result "vl" "Velocity"', f20.10, 2x, 'Vector OnGaussPoints "Cell Center"')
  120 FORMAT('Result "Regime" "Regime"', f20.10, 2x, 'Scalar OnGaussPoints "Cell Center"')   
!
      END SUBROUTINE GID_out_cell_lsj
