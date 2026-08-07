!
      SUBROUTINE GeoData_user_opr1000_1x1(nn)
!
!     Define user-defined geometry data for OPR1000 Full vessel (rod-scale)
!
      USE Zcore       , ONLY: np,myrank
      USE Znum_cell   , ONLY: i_neigh_tmp,j_neigh_tmp, &
                              sv_tmp1
      USE Zporous , ONLY: chn_type_tmp
      USE Zrv_ncell   , ONLY: asm_nx,asm_ny,asm_nz,asm_ni,asm_ni2,chn_nx,chn_ny
      USE Zrv_ncell   , ONLY: master_to_assem0,master_to_assem1

      USE Zrv_hts_2d  , ONLY:dz_fuel00
      USE MASTER4     , ONLY:NXY_TH,NZ_TH,NXYF,NPINX,nchn,num_asm,mst_to_asmi
      USE Zio_unit    , ONLY: unit_log
!
      IMPLICIT NONE

      !input
      INTEGER nn
      
      INTEGER i,ii,j
      INTEGER i0,z0,i1,i2,i3,i4,i5,i6,i7,i8
      INTEGER nz0
      INTEGER scin      
      
      ! MASTER mapping
      INTEGER ny
      INTEGER m,n,k,kk
      INTEGER iz,iy,is,ic
      INTEGER,ALLOCATABLE::ja_f(:),ia_r(:)
      ! assembly mapping
      INTEGER id0,ir0,m1,ic0,icupid
      
      !total table
      INTEGER in
      INTEGER,ALLOCATABLE::asm0_to_cell(:,:)
      !invrt table
      INTEGER,ALLOCATABLE::asmi_to_cell(:),asmz_to_cell(:)
!
      INTEGER :: itmp(9)
 
       DATA scin/1/
     
!   
!======================================================================
! READ OPR1000 GeoData
!======================================================================

      IF(myrank.eq.0) THEN
         OPEN(unit=1550,file='subchannel_info.in',status='old', iostat=scin)
      ENDIF ! myrank
      IF(np.gt.1) CALL broadcast_i1(scin)
      IF(scin.ne.0)then
         IF(myrank.eq.0) WRITE(* ,*)'Program was terminated due to lack of <subchannel_info.in>!'
         IF(myrank.eq.0) WRITE(unit_log,*)'Program was terminated due to lack of <subchannel_info.in>!'
         CALL finalize_mpi
         STOP
      ENDIF

      ALLOCATE(asm_nx(nn),asm_ny(nn),asm_nz(nn),asm_ni(nn),asm_ni2(nn))
      ALLOCATE(chn_nx(nn),chn_ny(nn))
      ALLOCATE(chn_type_tmp(nn))

      IF(myrank.eq.0) READ(1550,*)i0,z0
      IF(np.gt.1) CALL broadcast_i1(i0)
      IF(np.gt.1) CALL broadcast_i1(z0)
      nchn=int(dsqrt(dble(i0/z0/177)))
      
      !total table
      ALLOCATE(asm0_to_cell(177*nchn,z0))
      asm0_to_cell=0
      !invrt table
      ALLOCATE(asmi_to_cell(177*nchn))
      ALLOCATE(asmz_to_cell(177*nchn))
      asmi_to_cell=0
      asmz_to_cell=0
      
      DO i=1,nn
         asm_nx(i)=0
         asm_ny(i)=0
         asm_nz(i)=0
         asm_ni(i)=0
         asm_ni2(i)=0
         chn_nx(i)=0
         chn_ny(i)=0
         chn_type_tmp(i)=0
      ENDDO
      
      num_asm=0
      DO ii=1,i0
         IF(myrank.eq.0) THEN
            READ(1550,49)i,i1,i2,i3,i4,i5,i6,i7,i8
            IF(np.gt.1) THEN
               itmp(1)=i
               itmp(2)=i1
               itmp(3)=i2
               itmp(4)=i3
               itmp(5)=i4
               itmp(6)=i5
               itmp(7)=i6
               itmp(8)=i7
               itmp(9)=i8
            ENDIF
         ENDIF
         IF(np.gt.1) THEN
            CALL broadcast_i(itmp,9)
            i =itmp(1)
            i1=itmp(2)
            i2=itmp(3)
            i3=itmp(4)
            i4=itmp(5)
            i5=itmp(6)
            i6=itmp(7)
            i7=itmp(8)
            i8=itmp(9)
         ENDIF
         asm_nx(i)=i1
         asm_ny(i)=i2
         asm_nz(i)=i3
         asm_ni(i)=i4
         asm_ni2(i)=i5
         chn_nx(i)=i6
         chn_ny(i)=i7
         chn_type_tmp(i)=i8
         
         !total table
         asm0_to_cell(i4,i3)=i
         !invrt table
         IF(i3.eq.1)then
            asmi_to_cell(i4)=i
         ENDIF   
         !number of asm, num_asm
         IF(i4.gt.num_asm)num_asm=i4
      ENDDO
49 format(20(i8,1x))

     !nz0=z0-2
      nz0=z0
      ALLOCATE(dz_fuel00(nz0))
      IF(myrank.eq.0) THEN
         DO j=1,nz0
            READ(1550,*)dz_fuel00(j)
         ENDDO
      ENDIF
      IF(np.gt.1) CALL broadcast_r(dz_fuel00,nz0)

      !MASTER variables => Fixed size
      nxyf  =i0/z0
      nxy_th=(nxyf+64)*4
      nz_th =nz0
      npinx=16

!======================================================================         
!........Connectivity for MASTER
!======================================================================         
!
!.....Radial Mapping file
!     2x2assembly indexing
      ny=30
      ALLOCATE(ja_f(ny),ia_r(ny))
      ja_f=(/10,10, 18,18, 22,22, 26,26, 26,26, &   !5
             30,30, 30,30, 30,30, 30,30, 30,30, &   !10
             26,26, 26,26, 22,22, 18,18, 10,10   /) !15
      ia_r=(/34,12, 10, 8,  8, 8,  6, 4,  6, 8, &   !5
              6, 4,  4, 4,  4, 4,  4, 4,  4, 4, &   !10
              6, 8,  6, 4,  6, 8,  8, 8, 10,12   /) !15
!
!.....MASTER mapping, Assembly connectivity
      ALLOCATE(master_to_assem0(NXY_TH,NZ_TH)) ! master_to_channel
      ALLOCATE(master_to_assem1(NXY_TH,NZ_TH)) ! master_to_cupid
      master_to_assem0=0
      master_to_assem1=0
      
!      ALLOCATE(master_to_assem1_cell(NXY_TH,NZ_TH))
!      ALLOCATE(master_to_assem1_rod(NXY_TH,NZ_TH))
!      master_to_assem1_cell=0
!      master_to_assem1_rod=0
      
      ALLOCATE(mst_to_asmi(NXY_TH))
      mst_to_asmi=0
      
      ! ic --> channel numbering, nf_input(ic,nz0_2d)
      n=0
      ic0=0
      DO iy=1,30
         ir0=ia_r(iy)
         id0=ja_f(iy)
         is =n+ir0
         IF(iy.eq.1)is=ir0
         DO j=1,id0
            m=is+j
            IF(MOD(j,2).eq.0)then
               m1=INT(j/2)
            ELSE
               m1=INT(j/2) 
               m1=m1+1
            ENDIF
            ic=ic0+m1
            master_to_assem0(m,1)=ic
            mst_to_asmi(m)=ic
            master_to_assem1(m,1)=asmi_to_cell(ic)
         ENDDO
         IF(MOD(iy,2).eq.0)ic0=ic
         n=m
      ENDDO
      DEALLOCATE(ia_r,ja_f)
      
      ! Axial numbering
      DO i=1,NXY_TH
         DO k=1,NZ_TH-1
            i0=master_to_assem1(i,k)
            IF(i0.ne.0)then
               in=asm_ni(i0)
               iz=asm_nz(i0)
               i1=asm0_to_cell(in,iz+1)
               master_to_assem1(i,iz+1)=i1
            ENDIF   
         ENDDO
      ENDDO
             
      ! additional cupid cell numbering for REAL reflector (arbitrary cell is OK)
      i=asm0_to_cell(1,1)
      IF(myrank.eq.0) THEN
         DO j=i_neigh_tmp(i),i_neigh_tmp(i+1)-1
            IF(sv_tmp1(j,3).lt.0.d0)kk=j_neigh_tmp(j)
         ENDDO
      ENDIF
      IF(np.gt.1) CALL broadcast_i1(kk)
!      
      DO i=1,NXY_TH
         DO k=1,NZ_TH
            icupid=master_to_assem1(i,k)
            IF(icupid.eq.0)then
               master_to_assem1(i,k)=kk
            ENDIF
         ENDDO
      ENDDO

      END SUBROUTINE GeoData_user_opr1000_1x1
