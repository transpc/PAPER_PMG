      SUBROUTINE Prep_MG_GarL
! 
! * * * * * * * * * * * * * * * * * * * * * * * *  * * * * * * * * * !
! ---
      USE MD_geometry, ONLY: nnode
      USE MD_MG_coord, ONLY: ialv,nnods,coordc
      USE MD_parameter, ONLY: ndim,ndom
      USE MD_MPI_MG, ONLY: iintf,inodegl,nintfs,inbdc,nnsend_m,nnrecv_m,             &
                           ibdomc,isptc,irptc,isintfc,irintfc,isiac,iriac,           &
                           nnbds,spts,rpts,sintfs,rintfs,nbdoms
      USE MD_MG_matrix, ONLY: ia2,ja2,iac,jac,ias,jas,nnzs,nnzr
      USE MD_MG_index, ONLY: nlevel,n_GC
      USE MD_MPI, ONLY: myrank 
      USE MD_MG_Global_C , ONLY: i_dir,nlv_glo,nnodeG,nnodeGC,Ainv,nnzG,nnzGC,aluG, igather

!
!-------------------
      IMPLICIT NONE
	  
      INTEGER(4):: ilv,i,j,nd,nnzss,nnzrr
      INTEGER(4):: nnzt,nnz2,ncolf1,ncolc1,ntmp,nnode2,nintf2,nnode2gl
      INTEGER(4):: nnbd,nnsend,nnrecv,nmaxgl
      INTEGER(4),DIMENSION(:),ALLOCATABLE::nbdom,spt,rpt,sintf,rintf,sia,ria
!      INTEGER(4)::nbdom(ndom),spt(ndom),rpt(ndom),sintf(nnsend_m),rintf(nnrecv_m)
!      INTEGER(4)::sia(nnsend_m+1),ria(nnrecv_m+1)
      REAL(8),DIMENSION(:,:),ALLOCATABLE::coord2
      
      ALLOCATE(nbdom(ndom),spt(ndom),rpt(ndom),sintf(nnsend_m),rintf(nnrecv_m))
      ALLOCATE(sia(nnsend_m+1),ria(nnrecv_m+1))
! 
      nbdom = 0
      spt = 0
      rpt = 0
      sintf = 0
      rintf = 0
      sia = 0
      ria = 0
! --- 
!   
     IF(nlevel.EQ.2) GOTO 300
     
     DO ilv = 2,nlevel-1
      
! 
! step 1: preprocessing:
      
      nnbd = inbdc(ilv)  
      nbdom(1:nnbd) = ibdomc(1:nnbd,ilv)
      spt(1:nnbd+1) = isptc(1:nnbd+1,ilv)
      rpt(1:nnbd+1) = irptc(1:nnbd+1,ilv)

      nnsend = spt(nnbd+1)-1
      nnrecv = rpt(nnbd+1)-1
      nmaxgl = 1
      
! 
      nnsend = MAX(0,nnsend)
      nnrecv = MAX(0,nnrecv)
      
! 
   IF(nnsend.GT.nnsend_m.OR.nnrecv.GT.nnrecv_m) THEN
       WRITE(*,*)'PMG error: nnsend_m is small-2',nnsend,nnrecv,nnsend_m,nnrecv_m
       STOP
   ENDIF
      
      sintf(1:nnsend) = isintfc(1:nnsend,ilv)
      rintf(1:nnrecv) = irintfc(1:nnrecv,ilv)
!
      
! ia and ja:
      nintf2 = iintf(ilv)
      nnode2 = ialv(ilv+1)-ialv(ilv)
      nnode2gl = inodegl(ilv)
      
      ncolf1 = ialv(ilv+1) - ialv(1)
      ncolc1 = ncolf1 -nnode

      ntmp = ncolc1 - nnode2
      nnzt = iac(ntmp+1)-1
      
      nnz2 = iac(ncolc1+1)-iac(ntmp+1)
      
      ALLOCATE(ia2(nnode2+1),ja2(nnz2))
      
      ia2(1:nnode2+1)=iac(ntmp+1:ncolc1+1) - nnzt
      ja2 (1:nnz2) = jac(nnzt+1:nnzt+nnz2) - ntmp
      
! - - - - 
      sia(1) = 1
      Do i=1,nnsend
          j=sintf(i)
          nd = ia2(j+1)-ia2(j)
          sia(i+1)=sia(i)+nd
          nmaxgl = MAX(nmaxgl,nd)
      End do
      
      isiac(1:nnsend+1,ilv) = sia(1:nnsend+1)

      nnzss = sia(nnsend+1)-1
!     
      ria(1) = 1
      Do i=1,nnrecv
          j=rintf(i)
          nd = ia2(j+1)-ia2(j)
          ria(i+1)=ria(i)+nd
          nmaxgl = MAX(nmaxgl,nd)
      End do      
      
      iriac(1:nnrecv+1,ilv) = ria(1:nnrecv+1)
! check:
      IF((ria(nnrecv+1)-ria(1)).NE.(ia2(nnode2+1)-ia2(nintf2+1))) THEN
          WRITE(*,*)'PMG error in ria',myrank,ria(nnrecv+1)-ria(1),ia2(nnode2+1)-ia2(nintf2+1)
          STOP
      ENDIF      
      
      nnzrr = ria(nnrecv+1)-1
  !     
! step 2: adding the csr for external nodes:
! coord
      ALLOCATE(coord2(ndim,nnode2gl))
      i = SUM(inodegl(2:ilv-1))
      coord2(1:ndim,1:nnode2gl) = coordc(1:ndim,i+1:i+nnode2gl)
        
    IF(nnbd.NE.0) THEN
      CALL send_receive_csrc(ndim,nmaxgl,nnbd,nnode2,nnode2gl,spt,rpt,sintf,rintf,                &
              nnzss,nnzrr,nnsend,nnrecv,sia,ria,nbdom,coord2,nnz2,ia2,ja2,ndom,nnsend_m,nnrecv_m)	 
    ENDIF
! update ja:
      jac(nnzt+1:nnzt+nnz2) =  ja2 (1:nnz2) + ntmp   
      
      DEALLOCATE(ia2,ja2,coord2)
      
     ENDDO
     

300  CONTINUE

  IF(n_GC.EQ.1) THEN

  CALL imapGZ_coarse(nintfs,nnods,nnzs,ias,jas)

  ! for GC
  IF(nlv_glo.EQ.0) THEN
!   for direct solver
    IF(i_dir.NE.0) THEN
    ALLOCATE(Ainv(nnodeG,nnodeG))
    Ainv = 0.d0
    ELSE
        ALLOCATE(aluG(nnzG))
    ENDIF
    
  ELSE
    CALL PREP_GMG_global_coarse
    nnode2 = nnodeGC(nlv_glo)
    
    IF(i_dir.NE.0) THEN
    ALLOCATE(Ainv(nnode2,nnode2))
    Ainv = 0.d0      
    ELSE
        ALLOCATE(aluG(nnzGC(nlv_glo)))
    ENDIF
    
  ENDIF
!
  ELSE 
! for CG-Parallel - - - - - - - - - - - 
      ALLOCATE(nbdoms(ndom),spts(ndom),rpts(ndom),sintfs(nnsend_m),rintfs(nnrecv_m))
      ilv = nlevel
      nnbds = inbdc(ilv)  
      nbdoms(1:nnbds) = ibdomc(1:nnbds,ilv)
      spts(1:nnbds+1) = isptc(1:nnbds+1,ilv)
      rpts(1:nnbds+1) = irptc(1:nnbd+1,ilv)

      nnsend = spts(nnbds+1)-1
      nnrecv = rpts(nnbds+1)-1
! 
      nnsend = MAX(0,nnsend)
      nnrecv = MAX(0,nnrecv)
! 
      sintfs(1:nnsend) = isintfc(1:nnsend,ilv)
      rintfs(1:nnrecv) = irintfc(1:nnrecv,ilv)
! - - - - - - - - - - - - - - - - - - - - - - - - 
  ENDIF
  
! 
! for gather
     IF(igather.EQ.1) THEN
     CALL imap_GATHER(ndom,myrank,nintfs,nnzs,ias)
     ENDIF

! 
  DEALLOCATE(nbdom,spt,rpt,sintf,rintf,sia,ria)
  DEALLOCATE(coordc)
  
      RETURN
      END
      
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
