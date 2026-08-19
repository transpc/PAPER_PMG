      subroutine stiffness_MG
! ---
      USE md_geometry, only: nnode,nnodegl
      USE MD_matrix, only: nnz, ia, ja,ju, au
      USE MD_MG_matrix, ONLY: ia1,ja1, ia2, ja2, au1, au2, iac,jac,juc,auc,diagrc,        &
                              iar, jar, Xrest, iar2, jar2, Xrest2,                        &
                              nnzs,ias,jas,aus,nnz1,nnz2
      USE MD_MG_index, ONLY: nlevel,mxnbne,n_GC,il1_gs
	  USE MD_MG_coord, only: nnode1,nnode2,nnode1gl,nnode2gl,ialv,nnods,inmax
      USE MD_MPI, ONLY: nnbd,spt,rpt,sintf,rintf,nbdom,nnsend,nnrecv,nintf,myrank
      USE MD_MPI_MG, ONLY:  siaf,riaf,iintf,inodegl,nintf1,nintf2,nintfs
      use omp_lib	  
      
      implicit none  
!     
      integer(4) jmax,i,j,k,ilv,nnzt,ncolf1,ncolc1,ntmp,nnzr2,nnzt1,ntmpf
      REAL(8) tmp
      real*8 time_begin,time_end,time_cpu
      real*8 time_begin0,time_end0,time_cpu0
	  
!      jmax = mxnbne
     
 
!      call system_clock(i,j,k)
!      time_begin0=real(i,kind=8)/real(j,kind=8)  
	  
! initial fine level
      nnode1 = nnode
      nnz1 = nnz
      nnode1gl = nnodegl
      nintf1 = nintf
!      ALLOCATE(ia1(nnode1+1),ja1(nnz1),au1(nnz1))
!      ia1 = ia
!      ja1 = ja
!      au1 = au
      ALLOCATE(ia1(1),ja1(1),au1(1))
      ia1 = 0
      ja1 = 0
      au1 = 0.d0
     
! test-MG
!     write(*,*)'Gar-pre'
 
      DO ilv = 1, nlevel-1
          
! coarse level
      nintf2 = iintf(ilv+1)
      nnode2gl = inodegl(ilv+1)

      nnode2 = ialv(ilv+2)-ialv(ilv+1)
      
      ncolf1 = ialv(ilv+2) - ialv(1)
      ncolc1 = ncolf1 -nnode

      ntmp = ncolc1 - nnode2
      nnzt = iac(ntmp+1)-1
      
      nnz2 = iac(ncolc1+1)-iac(ntmp+1)
      
      ALLOCATE(ia2(nnode2+1),ja2(nnz2),au2(nnz2))
      
    !$omp PARALLEL
	
    !$omp DO
	  DO i=1, nnode2+1
      ia2(i)=iac(i+ntmp) - nnzt
	  ENDDO
    !$omp END DO
	  
    !$omp DO
	  DO i=1, nnz2
      ja2 (i) = jac(i+nnzt) - ntmp
	  ENDDO
    !$omp END DO
	
    !$omp END PARALLEL
      
! for Xrest:
      nnzt1 = iar(ntmp+1)-1
      nnzr2 = iar(ncolc1+1)-iar(ntmp+1)
      
      ntmpf = ialv(ilv+1) - ialv(1)-nnode1
      
      ALLOCATE(iar2(nnode2+1),jar2(nnzr2),Xrest2(nnzr2))

    !$omp PARALLEL
	
    !$omp DO
	  DO i=1, nnode2+1
       iar2(i)=iar(i+ntmp)- nnzt1
	  ENDDO
    !$omp END DO
	
    !$omp DO
	  DO i=1, nnzr2
       jar2 (i)=jar(i+nnzt1) - ntmpf
       Xrest2 (i)=Xrest(i+nnzt1)
	  ENDDO
    !$omp END DO
	
    !$omp END PARALLEL
!
       jmax = 4*inmax(ilv)
       
!! step 1: transfer data of Af for external row: 
     
    IF(nnbd.NE.0) THEN
        
      IF(ilv.EQ.1) THEN
       nnsend = spt(nnbd+1)-1
       nnrecv = rpt(nnbd+1)-1

      CALL send_receive_mtf(nnbd,nnsend,nnrecv,spt,rpt,sintf,rintf,          &
              siaf,riaf,nbdom,ia,au)
      ELSE
    
      CALL MD_S_R_MT(ilv,ia1,au1)

      ENDIF
      
    ENDIF
    

! step 2: Galerkin condition
      
! 
    IF(ilv.EQ.1) THEN
      CALL stiff_coarse_P(nnode1gl,nnode1,nintf1,nnode2,nnode2gl,nintf2,jmax,      &
            nnz1,nnz2,nnzr2,ia,ja,ia2,ja2,iar2,jar2,au,au2,Xrest2)
    ELSE
      CALL stiff_coarse_P(nnode1gl,nnode1,nintf1,nnode2,nnode2gl,nintf2,jmax,      &
            nnz1,nnz2,nnzr2,ia1,ja1,ia2,ja2,iar2,jar2,au1,au2,Xrest2)
    ENDIF
!
      DEALLOCATE(iar2,jar2,Xrest2)
       
! adding to array
    !$omp PARALLEL DO
	   DO i=1, nnz2
       auc (i+nnzt) = au2 (i)
	   ENDDO
    !$omp END PARALLEL DO
! update for fine level
       DEALLOCATE(ia1,ja1,au1)
       nnode1 = nnode2
       nnz1 = nnz2
       nintf1 = nintf2
       nnode1gl = nnode2gl
       
       ALLOCATE(ia1(nnode1+1),ja1(nnz1),au1(nnz1))

    !$omp PARALLEL
	
    !$omp DO
	  DO i=1, nnode1+1
       ia1(i) = ia2(i)
	  ENDDO
    !$omp END DO	   

    !$omp DO
	  DO i=1, nnz1	
       ja1(i) = ja2(i)
       au1(i) = au2(i)
	  ENDDO
    !$omp END DO

    !$omp END PARALLEL

       DEALLOCATE(ia2,ja2,au2)

!     FOr diagrc

        ntmp = ncolc1 - nnode2
    !$omp PARALLEL DO private(j,k,tmp)

        DO i = 1, nintf2
                 j = juc(ntmp+i)
            tmp = auc(j)
! il1_gs=1: 랭크 밖(고스트, 레벨 내 로컬번호 > nintf2) 열의 |a_ij| 를 대각에
! 가산 — l1-hybrid-GS (해결책 A). 고스트 열은 jac > ntmp+nintf2 로 판별
            IF (il1_gs .EQ. 1) THEN
               DO k = iac(ntmp+i), iac(ntmp+i+1)-1
                  IF (jac(k) .GT. ntmp+nintf2) tmp = tmp + ABS(auc(k))
               ENDDO
            ENDIF
            diagrc(ntmp+i) = 1.d0/tmp
        ENDDO

    !$omp END PARALLEL DO
        
      ENDDO
      
!
    
! for the coarsest level:
      IF(nnode1.NE.nnods) THEN
          WRITE(999,*)'error in nnods-stiffness-MG'
          STOP
      ENDIF
      
      
    !$omp PARALLEL DO
	
	  DO i=1, nnz1
      aus(i) = au1(i)
!      alus(i) = au1(i)
	  ENDDO

    !$omp END PARALLEL DO
!      CALL pc_ilu1(nnode1,i,ias,jas,jus,alus)
      
      DEALLOCATE(ia1,ja1,au1)
      
!
!      call system_clock(i,j,k)
!      time_begin=real(i,kind=8)/real(j,kind=8)  
      
      
      IF(n_GC.EQ.1) THEN
        CALL stiffness_GC_all(nintfs,nnods,nnzs,ias,aus)
        
      ENDIF
! 
      
!     IF(myrank.EQ.0) THEN

!      call system_clock(i,j,k)
!      time_end0=real(i,kind=8)/real(j,kind=8)
!          time_cpu0 = time_end0-time_begin0
!      write(101,*)'cpu Gar in=',time_cpu0

!          time_cpu = time_end0-time_begin
!      write(*,*)'cpu Gar-GC=',time_cpu
!     ENDIF
    
!     write(*,*)'Gar-post'
      return	  
    End
    
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
! Galerkin Formula: 
    
SUBROUTINE stiff_coarse_P(nnodegl,nnode,nintf,nnode1,nnode1gl,nintf1,nfmax,nnz, &
                              nnz1,nnzr,ia,ja,ia1,ja1,iar,jar,au,au1,Xr)
    
! ************************************************************************************!
! this subroutine calculates the coarse-stiff matrix by Galerkin formular             !
!     Ac = R*Af*R(T)                                                                  !
!   or: Ac(I,J) = R(I,k)*Af(k,l)*R(J,l)                                               !
!   for each I, J on coarse grid, we only find k,l on fine-grid such that:            !
!   R(I,k) =/ 0 and R(J,l) =/0                                                        !
! inlet: Af, R                                                                        !
! outlet: Ac                                                                          !
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! for each node I on coarse grid                                                      !
! 1- VI(1:nf) = R(I,k)*Af(k,l) : this is just matrix-vector(dense) multiplication     !
!     several nodes k,l                                                               !  
!  for all neibough coarse-node J of node I                                           !              
! 2- A(I,J) = VI(l)*R(J,l) : this is just vector-vector (dense) multiplication        !
!     several node l for each node J                                                  !
! ------------------------------------------------------------------------------------!
      
!----notes that: nnode-> fine, nnode1-> coarse
!----nfmax: maximum fine node in neighbor elements of coarse node I
      use omp_lib
      
      implicit none
! ---inlet
      integer nnode,nnode1,nfmax,nintf,nintf1,nnodegl,nnode1gl
      integer nnz,nnz1,nnzr
      integer ia(nnode+1),ja(nnz)
      integer ia1(nnode1+1),ja1(nnz1)
      integer iar(nnode1+1),jar(nnzr)
      real*8 au(nnz),Xr(nnzr)
!      
! ---outlet 
      real*8 au1(nnz1)
! ---temp
      integer i,j,k,l,i1,i2,j1,j2,id,imax,jmax,ll
      real*8 s
      INTEGER(4) ni(nfmax),nj(nfmax)
      REAL(8) pi(nfmax),pj(nfmax)
      real*8,dimension(:),allocatable::vi
      
! ----------------------------------------------------------!
      allocate(vi(nnodegl))
      
!      au1 = 0.d0
!      ni = 0
!      vi = 0.d0 
!      pi = 0.d0
!      nj = 0
!      pj = 0.d0
    !$omp PARALLEL
    !$omp DO
      DO i=1,nnz1 	
      au1(i) = 0.d0
	  ENDDO
	!$omp END DO
	
    !$omp DO
      DO i=1,nnodegl
      vi(i) = 0.d0 
      ENDDO
	!$omp END DO	  

    !$omp END PARALLEL 
! ---
    !$omp PARALLEL DO private(i1,i2,imax,ni,pi,nj,pj,id,k,l,ll,j,j1,j2,jmax,s) firstprivate(vi)
  
      DO i = 1, nintf1
! ---for node i              
         i1 = iar(i)
         i2 = iar(i+1)-1		 
         imax = i2-i1+1
         ni(1:imax) = jar(i1:i2)
         pi(1:imax) = Xr(i1:i2)    
!    
         do id = 1,imax
            k=ni(id)
            
            do l = ia(k),ia(k+1)-1
            ll=ja(l)
            
            vi(ll) = vi(ll) + pi(id)*au(l)

            end do
         end do		
         
! ---for node j
!        
         do id = ia1(i),ia1(i+1)-1
		    j = ja1(id)

            j1 = iar(j)
            j2 = iar(j+1)-1		 
            jmax = j2-j1+1
            nj(1:jmax) = jar(j1:j2)
            pj(1:jmax) = Xr(j1:j2)
!---
            s = 0.d0
            do l = 1,jmax
                 IF(nj(l).LE.0) then

                     cycle
                 endif
                
            s = s + vi(nj(l))*pj(l)
            end do
            au1(id) = s
 
         end do
! reset vi=0
         do id = 1,imax
             k=ni(id)

            do l = ia(k),ia(k+1)-1
            vi(ja(l)) = 0.0
            end do
         end do	
!
      END DO
      
    !$omp end PARALLEL DO 
      
! -------------------------------------------------------------!
      DEALLOCATE(vi)
      
      return
      
    END

! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !