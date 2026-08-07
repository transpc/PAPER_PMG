!
      SUBROUTINE udfn_subchannel_rod
!
!     Define USEr-defined geometry data for fluidic device
!
      USE Zmpi           , ONLY: jperm
      USE Zcore          , ONLY: myrank
      USE Ztimecon       , ONLY: time,itim
      USE Zparam         , ONLY: pi
      USE Zconst1        , ONLY: vv_prob
      USE Zrv_ncell      , ONLY: master_to_rod,ncell_fuel_rod,cupid_cell_hts2d, &
                                 assem_nx0,assem_ny0,assem_nz0, &
                                 dnbr_ce1
      USE Zrv_hts_2d     , ONLY: t_fuel,nz0_2d,dz_fuel0
      USE Zrv_hts_2d     , ONLY: nr_2d 
      USE Zrv_mpi        , ONLY: jperm_fuel_rod
      USE MASTER4        , ONLY: NXYF,NPINX,NZ_TH
      USE Zcoord2        , ONLY: cell_leng
      USE Zrv_model      , ONLY: rv_model
      USE Zrv_ncell      , ONLY: asm_nx,asm_ny,asm_nz,chn_nx,chn_ny
      USE Zio_unit       , ONLY: unit_log
!
      IMPLICIT NONE
!.....Local variables
      INTEGER i,j,k,l,p,q,ncell,nrod,nassem
      INTEGER ncall
      CHARACTER *30 s_ncall,s_name,s_name1,s_name2,s_ncell,s_p,s_q,s_itim
      LOGICAL, SAVE::guide
      REAL(8) point_height,m
!.....Local arrays
      REAL(8),ALLOCATABLE :: rod_tsolid_FC(:,:,:,:,:),rod_tsolid_max_FC(:,:,:,:,:)
      REAL(8),ALLOCATABLE :: cell_height(:)
      REAL(8),ALLOCATABLE :: rod_DNBR_CE_FC(:,:,:,:,:)
      !CUPID-RV
!.....Local variables
      INTEGER iz,ic,irod,icupid,icupid_j,ax,ay,az
      INTEGER cx,cy,ci,ji
!.....Local arrays
!     OPR1000 - ROD paraview
!.....Local variables
      REAL(8) rod_diameter,rod_pitch,rod_wall 
!
      DATA guide/.false./
!
! Variables
      ncall=0
      point_height=0.d0
      i=0
      j=0
      k=0
      l=0
!
!......Initial variables
!
      IF(vv_prob.eq.'single_assem')THEN
      ELSEIF(vv_prob.eq.'APR1400_fullcore')THEN
         !ncell=34      ! axial cell number
      ELSEIF(vv_prob.eq.'OPR1000_fullcore')THEN
      ELSEIF(vv_prob.eq.'OPR1000_quarter_core')THEN
      ELSEIF(vv_prob.eq.'OPR1000_fullcore_modmesh02_rv')THEN
         ncell=nz0_2d-1
         nrod=NPINX       ! 16x16 rods
         nassem=15
         guide=.true.
      ELSEIF(vv_prob.eq.'OPR1000_fullcore_modmesh02_rv_vessel')THEN
         ncell=NZ_TH-1
         nrod=NPINX       ! 16x16 rods
         nassem=15
         guide=.true.

         rod_diameter=0.00475d0
         rod_pitch=0.012852d0
         rod_wall=0.006426d0

      ENDIF


      ALLOCATE(rod_tsolid_FC(nrod*nrod,1,ncell,nassem,nassem))
      ALLOCATE(rod_tsolid_max_FC(nrod*nrod,1,ncell,nassem,nassem))
      ALLOCATE(rod_DNBR_CE_FC(nrod*nrod,1,ncell,nassem,nassem))
      ALLOCATE(cell_height(ncell))
      rod_tsolid_FC=0.0d0
      rod_tsolid_max_FC=0.0d0
      rod_DNBR_CE_FC=0.0d0
      cell_height=0.0d0


      ! CUPID-RV
      !Old mapping
if(0)then
      DO iz=1,ncell
         DO k=1,NXYF
            DO j=1,NPINX
               DO i=1,NPINX
                  ic=master_to_rod(i,j,k,iz)
                  DO irod=1,ncell_fuel_rod
                     m=jperm_fuel_rod(irod)
                     IF(ic.eq.m)then
                        icupid=cupid_cell_hts2d(irod)
                        icupid_j=jperm(icupid)
                        ax=assem_nx0(icupid_j)
                        ay=assem_ny0(icupid_j)
                        az=assem_nz0(icupid_j)
                        IF(rv_model.eq.1)then 
                           rod_tsolid_max_FC(NPINX*(j-1)+i,1,az,ax,ay)=t_fuel(irod,1)
                           rod_tsolid_FC(NPINX*(j-1)+i,1,az,ax,ay)    =t_fuel(irod,nr_2d)
                        ENDIF
                        rod_DNBR_CE_FC(NPINX*(j-1)+i,1,az,ax,ay)   =dnbr_ce1(irod)
 
                        cell_height(iz)=cell_leng(icupid,3)
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
endif

!new mapping
      DO irod=1,ncell_fuel_rod
         ic=cupid_cell_hts2d(irod)
         ji=jperm(ic)
         ax=asm_nx(ji)
         ay=asm_ny(ji)
         az=asm_nz(ji)
         cx=chn_nx(ji)
         cy=chn_ny(ji)
         ci=npinx*(cy-1)+cx
         rod_tsolid_max_FC(ci,1,az,ax,ay)=t_fuel(irod,1)
         rod_tsolid_FC(ci,1,az,ax,ay)    =t_fuel(irod,nr_2d)
         rod_DNBR_CE_FC(ci,1,az,ax,ay)   =dnbr_ce1(irod)
      ENDDO

!if(myrank.eq.0)then
!write(*,*)(cell_height(i),i=1,ncell)
!endif
!stop 'cell_height'



!      DO i=1,nassem
!         DO j=1,nassem
!            DO k=1,ncell
!               CALL allreduce_r(rod_tsolid_max_FC(:,1,k,j,i),dum1(:,1,k,j,i),NPINX*NPINX)
!               CALL allreduce_r(rod_tsolid_FC(:,1,k,j,i),dum1(:,1,k,j,i),NPINX*NPINX)
!               CALL allreduce_r(rod_DNBR_CE_FC(:,1,k,j,i),dum3(:,1,k,j,i),NPINX*NPINX)
!            ENDDO
!         ENDDO
!      ENDDO
!      rod_tsolid_max_FC=dum1
!      rod_tsolid_FC=dum2
!      rod_DNBR_CE_FC=dum3
      

      IF(rv_model.eq.1)then
         CALL allreducei_r(rod_tsolid_max_FC(1,1,1,1,1),NPINX*NPINX*ncell*nassem*nassem)
         CALL allreducei_r(rod_tsolid_FC(1,1,1,1,1),NPINX*NPINX*ncell*nassem*nassem)
      ENDIF

      CALL allreducei_r(rod_DNBR_CE_FC(1,1,1,1,1),NPINX*NPINX*ncell*nassem*nassem)

!      CALL allreduce_r(cell_height,cell_height0,ncell)
      cell_height(1:ncell)=dz_fuel0(1:ncell)
!if(myrank.eq.0)then
!write(*,*)(cell_height(i),i=1,ncell)
!endif
!stop 'cell_height'

      IF(myrank.eq.0)write(*,*)'Allreduce rod_tsolid_max_FC/rod_tsolid_FC/rod_DNBR_CE_FC to myrank=0'



!
IF(1)THEN
!
!......OPR1000_fullcore
!
      IF(vv_prob.eq.'OPR1000_fullcore_modmesh02_rv' .or. &
         vv_prob.eq.'OPR1000_fullcore_modmesh02_rv_vessel')THEN
         IF(myrank.eq.0)THEN
            write(*,*)  "set vtk file at time=",time,itim
            WRITE(unit_log,*) "set vtk file at time=",time,itim
            ncall=ncall+1
            DO p=1,15
               DO q=1,15
                  IF(p.eq.1.and.q.le.5 .or. p.eq.1.and.q.ge.11 .or. p.eq.15.and.q.le.5 .or. p.eq.15.and.q.ge.11) cycle
                  IF(p.eq.2.and.q.le.3 .or. p.eq.2.and.q.ge.13 .or. p.eq.14.and.q.le.3 .or. p.eq.14.and.q.ge.13) cycle
                  IF(p.eq.3.and.q.le.2 .or. p.eq.3.and.q.ge.14 .or. p.eq.13.and.q.le.2 .or. p.eq.13.and.q.ge.14) cycle
                  IF(p.eq.4.and.q.le.1 .or. p.eq.4.and.q.ge.15 .or. p.eq.12.and.q.le.1 .or. p.eq.12.and.q.ge.15) cycle
                  IF(p.eq.5.and.q.le.1 .or. p.eq.5.and.q.ge.15 .or. p.eq.11.and.q.le.1 .or. p.eq.11.and.q.ge.15) cycle

                  WRITE( s_p,'(I0)') p
                  s_name1='x'
                  WRITE( s_q,'(I0)') q
                  s_name2='.vtk.'
                  WRITE( s_ncall,'(I0)') ncall
                  WRITE( s_itim,'(I0)') itim
                 !s_name='y.vtk/'//trim(s_p)//trim(s_name1)//trim(s_q)//trim(s_name2)//trim(s_ncall)
                  s_name='y.vtk/'//trim(s_p)//trim(s_name1)//trim(s_q)//trim(s_name2)//trim(s_itim)
                  
                  OPEN (3,file=s_name)
                  
                  write(3,'(a)')'# vtk DataFile Version 2.0'
                  write(3,'(a)')'Rod Data'
                  write(3,'(a)')'ASCII'
                  write(3,'(a)')'DATASET UNSTRUCTURED_GRID'
                  WRITE( s_ncell,'(I0)') 8*4*ncell*nrod*nrod
                  write(3,'(a)')'POINTS       '//trim(s_ncell)//' double'           ! final line number of points-5
                  
                  DO j=1,nrod         ! y-direction             1 2
                     DO i=1,nrod      ! x-direction             3 4
                        IF(guide)THEN
                        IF(i.eq.4.or.i.eq.5.or.i.eq.8.or.i.eq.9.or.i.eq.12.or.i.eq.13)THEN
                           IF(j.eq.4.or.j.eq.5.or.j.eq.8.or.j.eq.9.or.j.eq.12.or.j.eq.13)THEN
                              IF(MOD(i+j,8).le.2)THEN       ! guide tube
                                 DO k=1,2                    ! guide tube(1,:), guide tube(2,:)   : upper part
                                    point_height=0.d0
                                    DO l=1,ncell
                                       write(3,*)    (rod_wall-rod_wall)+rod_wall*(k-1)+rod_pitch*(i-1)+0.5d0*rod_wall*(1-mod(j,4))*(1-mod(i,4))*(2-k)+(p-1)*0.20773d0+0.00213d0, &
                                          (rod_wall+(nrod-1)*rod_pitch+rod_wall)-rod_pitch*(j-1)-0.25d0*rod_wall*(1-mod(j,4))-((-1)**mod(i,4))*0.25d0*rod_wall*(2-k)*(1-mod(j,4))+(15-q)*0.20773d0+0.00213d0, point_height
                                       write(3,*)     rod_wall          +rod_wall*(k-1)+rod_pitch*(i-1)-0.5d0*rod_wall*(1-mod(j,4))*(mod(i,4))*(k-1)+(p-1)*0.20773d0+0.00213d0, &
                                          (rod_wall+(nrod-1)*rod_pitch+rod_wall)-rod_pitch*(j-1)-0.25d0*rod_wall*(1-mod(j,4))+((-1)**mod(i,4))*0.25d0*rod_wall*(k-1)*(1-mod(j,4))+(15-q)*0.20773d0+0.00213d0, point_height
                                       write(3,*)    (rod_wall-rod_wall)+rod_wall*(k-1)+rod_pitch*(i-1)+(0.25d0*rod_wall-0.25d0*rod_wall*(k-1))*(1-mod(i,4))+(p-1)*0.20773d0+0.00213d0, &
                                           rod_wall+(nrod-1)*rod_pitch          -rod_pitch*(j-1)+(15-q)*0.20773d0+0.00213d0,point_height
                                       write(3,*)     rod_wall+rod_wall*(k-1)-0.25d0*rod_wall*(k-1)*mod(i,4)+rod_pitch*(i-1)+(p-1)*0.20773d0+0.00213d0, &
                                          rod_wall+(nrod-1)*rod_pitch           -rod_pitch*(j-1)+(15-q)*0.20773d0+0.00213d0,point_height

                                       point_height=point_height+cell_height(l) !180205
                                       !
                                       write(3,*)    (rod_wall-rod_wall)+rod_wall*(k-1)+rod_pitch*(i-1)+0.5d0*rod_wall*(1-mod(j,4))*(1-mod(i,4))*(2-k)+(p-1)*0.20773d0+0.00213d0, &
                                          (rod_wall+(nrod-1)*rod_pitch+rod_wall)-rod_pitch*(j-1)-0.25d0*rod_wall*(1-mod(j,4))-((-1)**mod(i,4))*0.25d0*rod_wall*(2-k)*(1-mod(j,4))+(15-q)*0.20773d0+0.00213d0, point_height
                                       write(3,*)     rod_wall          +rod_wall*(k-1)+rod_pitch*(i-1)-0.5d0*rod_wall*(1-mod(j,4))*(mod(i,4))*(k-1)+(p-1)*0.20773d0+0.00213d0, &
                                          (rod_wall+(nrod-1)*rod_pitch+rod_wall)-rod_pitch*(j-1)-0.25d0*rod_wall*(1-mod(j,4))+((-1)**mod(i,4))*0.25d0*rod_wall*(k-1)*(1-mod(j,4))+(15-q)*0.20773d0+0.00213d0, point_height
                                       write(3,*)    (rod_wall-rod_wall)+rod_wall*(k-1)+rod_pitch*(i-1)+(0.25d0*rod_wall-0.25d0*rod_wall*(k-1))*(1-mod(i,4))+(p-1)*0.20773d0+0.00213d0, &
                                           rod_wall+(nrod-1)*rod_pitch          -rod_pitch*(j-1)+(15-q)*0.20773d0+0.00213d0,point_height
                                       write(3,*)     rod_wall+rod_wall*(k-1)-0.25d0*rod_wall*(k-1)*mod(i,4)+rod_pitch*(i-1)+(p-1)*0.20773d0+0.00213d0, &
                                          rod_wall+(nrod-1)*rod_pitch           -rod_pitch*(j-1)+(15-q)*0.20773d0+0.00213d0,point_height
                                    ENDDO
                                 ENDDO
                                 DO k=2,1,-1                    ! guide tube(3,:), guide tube(4,:)
                                    point_height=0.d0
                                    DO l=1,ncell
                                       write(3,*)     (rod_wall-rod_wall)+rod_wall*(k-1)+rod_pitch*(i-1)+0.25d0*rod_wall*(1-mod(i,4))*(2-k)+(p-1)*0.20773d0+0.00213d0, &
                                          rod_wall+(nrod-1)*rod_pitch-rod_pitch*(j-1)+(15-q)*0.20773d0+0.00213d0,point_height       ! guide tube(k,1)
                                       write(3,*)     rod_wall           +rod_wall*(k-1)+rod_pitch*(i-1)-0.25d0*rod_wall*mod(i,4)*(k-1)+(p-1)*0.20773d0+0.00213d0, &
                                          rod_wall+(nrod-1)*rod_pitch-rod_pitch*(j-1)+(15-q)*0.20773d0+0.00213d0,point_height       ! guide tube(k,2)
                                       write(3,*)     (rod_wall-rod_wall)+rod_wall*(k-1)+rod_pitch*(i-1)+0.5d0*rod_wall*(1-mod(i,4))*mod(j,4)*(2-k)+(p-1)*0.20773d0+0.00213d0, &
                                          (rod_wall+(nrod-1)*rod_pitch-rod_wall)-rod_pitch*(j-1)+0.25d0*rod_wall*mod(j,4)+((-1)**mod(i,4))*0.25d0*rod_wall*(2-k)*mod(j,4)+(15-q)*0.20773d0+0.00213d0,point_height
                                       write(3,*)     rod_wall           +rod_wall*(k-1)+rod_pitch*(i-1)-0.5d0*rod_wall*mod(i,4)*mod(j,4)*(k-1)+(p-1)*0.20773d0+0.00213d00, &
                                          (rod_wall+(nrod-1)*rod_pitch-rod_wall)-rod_pitch*(j-1)+0.25d0*rod_wall*mod(j,4)-((-1)**mod(i,4))*0.25d0*rod_wall*(k-1)*mod(j,4)+(15-q)*0.20773d0+0.00213d0,point_height
                  
                                       point_height=point_height+cell_height(l) !180205
                                       !
                                       write(3,*)     (rod_wall-rod_wall)+rod_wall*(k-1)+rod_pitch*(i-1)+0.25d0*rod_wall*(1-mod(i,4))*(2-k)+(p-1)*0.20773d0+0.00213d0, &
                                          rod_wall+(nrod-1)*rod_pitch-rod_pitch*(j-1)+(15-q)*0.20773d0+0.00213d0,point_height       ! guide tube(k,1)
                                       write(3,*)     rod_wall           +rod_wall*(k-1)+rod_pitch*(i-1)-0.25d0*rod_wall*mod(i,4)*(k-1)+(p-1)*0.20773d0+0.00213d0, &
                                          rod_wall+(nrod-1)*rod_pitch-rod_pitch*(j-1)+(15-q)*0.20773d0+0.00213d0,point_height       ! guide tube(k,2)
                                       write(3,*)     (rod_wall-rod_wall)+rod_wall*(k-1)+rod_pitch*(i-1)+0.5d0*rod_wall*(1-mod(i,4))*mod(j,4)*(2-k)+(p-1)*0.20773d0+0.00213d0, &
                                          (rod_wall+(nrod-1)*rod_pitch-rod_wall)-rod_pitch*(j-1)+0.25d0*rod_wall*mod(j,4)+((-1)**mod(i,4))*0.25d0*rod_wall*(2-k)*mod(j,4)+(15-q)*0.20773d0+0.00213d0,point_height
                                       write(3,*)     rod_wall           +rod_wall*(k-1)+rod_pitch*(i-1)-0.5d0*rod_wall*mod(i,4)*mod(j,4)*(k-1)+(p-1)*0.20773d0+0.00213d0, &
                                          (rod_wall+(nrod-1)*rod_pitch-rod_wall)-rod_pitch*(j-1)+0.25d0*rod_wall*mod(j,4)-((-1)**mod(i,4))*0.25d0*rod_wall*(k-1)*mod(j,4)+(15-q)*0.20773d0+0.00213d0,point_height
                                    ENDDO
                                 ENDDO
                                 CYCLE
                              ENDIF
                           ENDIF
                        ENDIF
                        ENDIF
                  
                        DO k=1,2    ! normal rod
                           point_height=0.d0
                           DO l=1,ncell
                              write(3,*)     (rod_wall-0.75*rod_diameter) +0.75*rod_diameter*(k-1) +rod_pitch*(i-1)+(p-1)*0.20773d0+0.00213d0, &
                                 (rod_wall+(nrod-1)*rod_pitch+0.75*rod_diameter) +0.25*rod_diameter*(k-1)-rod_pitch*(j-1)+(15-q)*0.20773d0+0.00213d0,point_height
                              write(3,*)     rod_wall                     +0.75*rod_diameter*(k-1) +rod_pitch*(i-1)+(p-1)*0.20773d0+0.00213d0, &
                                 (rod_wall+(nrod-1)*rod_pitch+rod_diameter)      -0.25*rod_diameter*(k-1)-rod_pitch*(j-1)+(15-q)*0.20773d0+0.00213d0,point_height
                              write(3,*)     (rod_wall-rod_diameter)+rod_diameter*(k-1)+rod_pitch*(i-1)+(p-1)*0.20773d0+0.00213d0, &
                                 rod_wall+(nrod-1)*rod_pitch                                             -rod_pitch*(j-1)+(15-q)*0.20773d0+0.00213d0,point_height
                              write(3,*)     rod_wall            +rod_diameter*(k-1)+rod_pitch*(i-1)+(p-1)*0.20773d0+0.00213d0, &
                                 rod_wall+(nrod-1)*rod_pitch            -rod_pitch*(j-1)+(15-q)*0.20773d0+0.00213d0,point_height
                              !
                              point_height=point_height+cell_height(l) !180205
                              !
                              write(3,*)     (rod_wall-0.75*rod_diameter) +0.75*rod_diameter*(k-1) +rod_pitch*(i-1)+(p-1)*0.20773d0+0.00213d0, &
                                 (rod_wall+(nrod-1)*rod_pitch+0.75*rod_diameter) +0.25*rod_diameter*(k-1)-rod_pitch*(j-1)+(15-q)*0.20773d0+0.00213d0,point_height
                              write(3,*)     rod_wall                     +0.75*rod_diameter*(k-1) +rod_pitch*(i-1)+(p-1)*0.20773d0+0.00213d0, &
                                 (rod_wall+(nrod-1)*rod_pitch+rod_diameter)      -0.25*rod_diameter*(k-1)-rod_pitch*(j-1)+(15-q)*0.20773d0+0.00213d0,point_height
                              write(3,*)     (rod_wall-rod_diameter)+rod_diameter*(k-1)+rod_pitch*(i-1)+(p-1)*0.20773d0+0.00213d0, &
                                 rod_wall+(nrod-1)*rod_pitch                                             -rod_pitch*(j-1)+(15-q)*0.20773d0+0.00213d0,point_height
                              write(3,*)     rod_wall            +rod_diameter*(k-1)+rod_pitch*(i-1)+(p-1)*0.20773d0+0.00213d0, &
                                 rod_wall+(nrod-1)*rod_pitch            -rod_pitch*(j-1)+(15-q)*0.20773d0+0.00213d0,point_height
                           ENDDO
                        ENDDO
                        DO k=2,1,-1
                           point_height=0.d0
                           DO l=1,ncell
                              write(3,*)     (rod_wall-rod_diameter)+rod_diameter*(k-1)+rod_pitch*(i-1)+(p-1)*0.20773d0+0.00213d0, &
                                               rod_wall+(nrod-1)*rod_pitch            -rod_pitch*(j-1)+(15-q)*0.20773d0+0.00213d0,point_height
                              write(3,*)     rod_wall            +rod_diameter*(k-1)+rod_pitch*(i-1)+(p-1)*0.20773d0+0.00213d0, &
                                               rod_wall+(nrod-1)*rod_pitch            -rod_pitch*(j-1)+(15-q)*0.20773d0+0.00213d0,point_height
                              write(3,*)     (rod_wall-0.75*rod_diameter) +0.75*rod_diameter*(k-1) +rod_pitch*(i-1)+(p-1)*0.20773d0+0.00213d0, &
                                              (rod_wall+(nrod-1)*rod_pitch-0.75*rod_diameter) -0.25*rod_diameter*(k-1)-rod_pitch*(j-1)+(15-q)*0.20773d0+0.00213d0,point_height
                              write(3,*)     rod_wall            +0.75*rod_diameter*(k-1) +rod_pitch*(i-1)+(p-1)*0.20773d0+0.00213d0, &
                                              (rod_wall+(nrod-1)*rod_pitch-rod_diameter)+0.25*rod_diameter*(k-1)-rod_pitch*(j-1)+(15-q)*0.20773d0+0.00213d0,point_height

                              point_height=point_height+cell_height(l) !180205
                     
                              write(3,*)     (rod_wall-rod_diameter)+rod_diameter*(k-1)+rod_pitch*(i-1)+(p-1)*0.20773d0+0.00213d0, &
                                               rod_wall+(nrod-1)*rod_pitch            -rod_pitch*(j-1)+(15-q)*0.20773d0+0.00213d0,point_height
                              write(3,*)     rod_wall            +rod_diameter*(k-1)+rod_pitch*(i-1)+(p-1)*0.20773d0+0.00213d0, &
                                               rod_wall+(nrod-1)*rod_pitch            -rod_pitch*(j-1)+(15-q)*0.20773d0+0.00213d0,point_height
                              write(3,*)     (rod_wall-0.75*rod_diameter) +0.75*rod_diameter*(k-1) +rod_pitch*(i-1)+(p-1)*0.20773d0+0.00213d0, &
                                              (rod_wall+(nrod-1)*rod_pitch-0.75*rod_diameter) -0.25*rod_diameter*(k-1)-rod_pitch*(j-1)+(15-q)*0.20773d0+0.00213d0,point_height
                              write(3,*)     rod_wall            +0.75*rod_diameter*(k-1) +rod_pitch*(i-1)+(p-1)*0.20773d0+0.00213d0, &
                                              (rod_wall+(nrod-1)*rod_pitch-rod_diameter)+0.25*rod_diameter*(k-1)-rod_pitch*(j-1)+(15-q)*0.20773d0+0.00213d0,point_height
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
                  
                  write(3,*)    'CELLS',ncell*4*nrod*nrod,(8*4*ncell*nrod*nrod) + (ncell*4*nrod*nrod)
                  DO j=1,nrod
                     DO i=1,nrod
                        DO k=1,4
                           DO l=1,ncell
                              write(3,20)     8,0+8*(l-1)+8*ncell*(k-1)+8*ncell*4*(i-1)+8*ncell*4*nrod*(j-1),1+8*(l-1)+8*ncell*(k-1)+8*ncell*4*(i-1)+8*ncell*4*nrod*(j-1),2+8*(l-1)+8*ncell*(k-1)+8*ncell*4*(i-1)+8*ncell*4*nrod*(j-1),3+8*(l-1)+8*ncell*(k-1)+8*ncell*4*(i-1)+8*ncell*4*nrod*(j-1),4+8*(l-1)+8*ncell*(k-1)+8*ncell*4*(i-1)+8*ncell*4*nrod*(j-1),5+8*(l-1)+8*ncell*(k-1)+8*ncell*4*(i-1)+8*ncell*4*nrod*(j-1),6+8*(l-1)+8*ncell*(k-1)+8*ncell*4*(i-1)+8*ncell*4*nrod*(j-1),7+8*(l-1)+8*ncell*(k-1)+8*ncell*4*(i-1)+8*ncell*4*nrod*(j-1)
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
                  write(3,*)    'CELL_TYPES',ncell*4*nrod*nrod
                  DO j=1,nrod
                     DO i=1,nrod
                        DO k=1,4
                           DO l=1,ncell
                              write(3,*)    11
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
                  write(3,*)    'CELL_DATA',ncell*4*nrod*nrod

              IF(rv_model.eq.1)then
                  write(3,'(a)')    'SCALARS Cladding_surface_temperature float'
                  write(3,'(a)')    'LOOKUP_TABLE default'
                  

                  
              
                  IF(p.ge.8)THEN
                     IF(q.ge.8)THEN
                        DO j=1,nrod
                           DO i=1,nrod
                              DO k=1,4
                                 DO l=1,ncell
                                    write(3,*)     rod_tsolid_FC(nrod*(j-1)+i,1,l,p,q)
                                 ENDDO
                              ENDDO
                           ENDDO
                        ENDDO
                     ELSEIF(q.le.7)THEN
                        DO j=1,nrod
                           DO i=1,nrod
                              !DO k=1,4
                              DO k=4,1,-1
                                 DO l=1,ncell
!                                    write(3,*)     rod_tsolid_FC(nrod*(16-j)+i,1,l,p,q) !snu-origin
                                    write(3,*)     rod_tsolid_FC(nrod*(j-1)+i,1,l,p,q) !jrlee
                                 ENDDO
                              ENDDO
                           ENDDO
                        ENDDO
                     ENDIF
                  ELSEIF(p.le.7)THEN
                     IF(q.ge.8)THEN
                        DO j=1,nrod
                           DO i=1,nrod
                              !DO k=1,4
                              DO k=2,1,-1
                                 DO l=1,ncell
!                                    write(3,*)     rod_tsolid_FC(nrod*(j-1)+(17-i),1,l,p,q) !snu-origin
                                    write(3,*) rod_tsolid_FC(nrod*(j-1)+i,1,l,p,q) !jrlee
                                 ENDDO
                              ENDDO
                              DO k=4,3,-1
                                 DO l=1,ncell
!                                    write(3,*)     rod_tsolid_FC(nrod*(j-1)+(17-i),1,l,p,q) !snu-origin
                                    write(3,*) rod_tsolid_FC(nrod*(j-1)+i,1,l,p,q) !jrlee
                                 ENDDO
                              ENDDO
                           ENDDO
                        ENDDO
                     ELSEIF(q.le.7)THEN
                        DO j=1,nrod
                           DO i=1,nrod
                              !DO k=1,4
                              DO k=3,4
                                 DO l=1,ncell
!                                    write(3,*)     rod_tsolid_FC(nrod*(16-j)+(17-i),1,l,p,q) !snu-origin
                                    write(3,*) rod_tsolid_FC(nrod*(j-1)+i,1,l,p,q) !jrlee
                                 ENDDO
                              ENDDO
                              DO k=1,2
                                 DO l=1,ncell
!                                    write(3,*)     rod_tsolid_FC(nrod*(16-j)+(17-i),1,l,p,q) !snu-origin
                                    write(3,*) rod_tsolid_FC(nrod*(j-1)+i,1,l,p,q) !jrlee
                                 ENDDO
                              ENDDO
                           ENDDO
                        ENDDO
                     ENDIF
                  ENDIF

                  write(3,'(a)')    'SCALARS Pellet_centerline_temperature float'
                  write(3,'(a)')    'LOOKUP_TABLE default'
              
                  IF(p.ge.8)THEN
                     IF(q.ge.8)THEN
                        DO j=1,nrod
                           DO i=1,nrod
                              DO k=1,4
                                 DO l=1,ncell
                                    write(3,*)     rod_tsolid_max_FC(nrod*(j-1)+i,1,l,p,q)
                                 ENDDO
                              ENDDO
                           ENDDO
                        ENDDO
                     ELSEIF(q.le.7)THEN
                        DO j=1,nrod
                           DO i=1,nrod
                              !DO k=1,4
                              DO k=4,1,-1
                                 DO l=1,ncell
                                    write(3,*)     rod_tsolid_max_FC(nrod*(j-1)+i,1,l,p,q)
                                 ENDDO
                              ENDDO
                           ENDDO
                        ENDDO
                     ENDIF
                  ELSEIF(p.le.7)THEN
                     IF(q.ge.8)THEN
                        DO j=1,nrod
                           DO i=1,nrod
                              !DO k=1,4
                              DO k=2,1,-1
                                 DO l=1,ncell
                                    write(3,*)     rod_tsolid_max_FC(nrod*(j-1)+i,1,l,p,q)
                                 ENDDO
                              ENDDO
                              DO k=4,3,-1
                                 DO l=1,ncell
                                    write(3,*)     rod_tsolid_max_FC(nrod*(j-1)+i,1,l,p,q)
                                 ENDDO
                              ENDDO
                           ENDDO
                        ENDDO
                     ELSEIF(q.le.7)THEN
                        DO j=1,nrod
                           DO i=1,nrod
                              !DO k=1,4
                              DO k=3,4
                                 DO l=1,ncell
                                    write(3,*)     rod_tsolid_max_FC(nrod*(j-1)+i,1,l,p,q)
                                 ENDDO
                              ENDDO
                              DO k=1,2
                                 DO l=1,ncell
                                    write(3,*)     rod_tsolid_max_FC(nrod*(j-1)+i,1,l,p,q)
                                 ENDDO
                              ENDDO
                           ENDDO
                        ENDDO
                     ENDIF
                  ENDIF
              endif
!!!!DNBR_CE
                  write(3,'(a)')    'SCALARS DNBR_CE float'
                  write(3,'(a)')    'LOOKUP_TABLE default'
              
                  IF(p.ge.8)THEN
                     IF(q.ge.8)THEN
                        DO j=1,nrod
                           DO i=1,nrod
                              DO k=1,4
                                 DO l=1,ncell
                                    write(3,*)     rod_DNBR_CE_FC(nrod*(j-1)+i,1,l,p,q)
                                 ENDDO
                              ENDDO
                           ENDDO
                        ENDDO
                     ELSEIF(q.le.7)THEN
                        DO j=1,nrod
                           DO i=1,nrod
                              !DO k=1,4
                              DO k=4,1,-1
                                 DO l=1,ncell
                                    write(3,*)     rod_DNBR_CE_FC(nrod*(j-1)+i,1,l,p,q)
                                 ENDDO
                              ENDDO
                           ENDDO
                        ENDDO
                     ENDIF
                  ELSEIF(p.le.7)THEN
                     IF(q.ge.8)THEN
                        DO j=1,nrod
                           DO i=1,nrod
                              !DO k=1,4
                              DO k=2,1,-1
                                 DO l=1,ncell
                                    write(3,*)     rod_DNBR_CE_FC(nrod*(j-1)+i,1,l,p,q)
                                 ENDDO
                              ENDDO
                              DO k=4,3,-1
                                 DO l=1,ncell
                                    write(3,*)     rod_DNBR_CE_FC(nrod*(j-1)+i,1,l,p,q)
                                 ENDDO
                              ENDDO
                           ENDDO
                        ENDDO
                     ELSEIF(q.le.7)THEN
                        DO j=1,nrod
                           DO i=1,nrod
                              !DO k=1,4
                              DO k=3,4
                                 DO l=1,ncell
                                    write(3,*)     rod_DNBR_CE_FC(nrod*(j-1)+i,1,l,p,q)
                                 ENDDO
                              ENDDO
                              DO k=1,2
                                 DO l=1,ncell
                                    write(3,*)     rod_DNBR_CE_FC(nrod*(j-1)+i,1,l,p,q)
                                 ENDDO
                              ENDDO
                           ENDDO
                        ENDDO
                     ENDIF
                  ENDIF
                  CLOSE(3)
               ENDDO
            ENDDO
         ENDIF
      ENDIF
            DEALLOCATE(rod_tsolid_FC,rod_tsolid_max_FC,rod_DNBR_CE_FC,cell_height)

ENDIF

20    FORMAT(9i10)

      RETURN
      END SUBROUTINE udfn_subchannel_rod
