 SUBROUTINE stiffness_GC(nintf,nnode,nnz,ia,au)
          
      USE MD_MG_Global_C, ONLY: i_dir,nlv_glo,nnodeG,nnzG,eG,rG,rG0,imapG,iaG,jaG,juG,auG,auG0,imapGZ, &
                                igather, nsengatA, irevgatA, idispA , imapgatA
      USE MD_MPI, ONLY: myrank
      USE MD_parameter, ONLY: ndom,ipar
      
      IMPLICIT NONE
          
!DEC$IF defined (mpi_flag)
      INCLUDE 'mpif.h'
!DEC$ENDIF
!
! input: 
  INTEGER nintf,nnode,nnz
  INTEGER ia(nnode+1)
  REAL*8 au(nnz)
! out:
! temp:
  INTEGER i,j,ierr,i1,i2,k,nd,ip

  REAL*8 time_begin,time_end,time_cpu 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
!step-1: transfer from local to global
        auG = 0.d0
        auG0 = 0.d0
         
    IF(igather.EQ.0) THEN
        
         
         DO i = 1,nintf
            i1 = ia(i)
                i2 = ia(i+1)
                 
            DO j=i1,i2-1
                auG0(imapGZ(j)) = au(j)
                
            ENDDO
        
            
         ENDDO   
!

    ENDIF
    
!          
         IF(ndom.EQ.1) THEN
             
             IF(igather.EQ.0) THEN
             auG = auG0
             ELSE
             auG(imapgatA(1:nnzG)) = au(1:nnzG)
             ENDIF
             
         ELSE
             
!DEC$IF defined (mpi_flag)
!    call mpi_barrier(mpi_comm_world,ierr)

! step-2: S&R to all processors

!      call system_clock(i,j,k)
!      time_begin=real(i,kind=8)/real(j,kind=8) 

        IF(igather.EQ.0) THEN
!       CALL MPI_ALLREDUCE(auG0,auG,nnzG,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
        CALL MPI_REDUCE(auG0,auG,nnzG,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)

        ELSE
           
    CALL MPI_GATHERV(au,nsengatA,MPI_double_precision,auG0,irevgatA,idispA,mpi_double_precision,0,mpi_comm_world, ierr)
        
!    CALL MPI_ALLGATHER(ausm,nnzsm,mpi_double_precision,auGAT,nnzsm,mpi_double_precision,mpi_comm_world, ierr)

        IF(myrank.EQ.0) THEN 
!         nd = 0
!          DO i = 1,nnzG
            
            DO j = 1,nnzG  
            auG(imapgatA(j)) = auG0(j)
            ENDDO
!          nd = nd+innzsm(ip)

!           ENDDO
         ENDIF

            
      ENDIF

                    
!    call mpi_barrier(mpi_comm_world,ierr)

!     IF(myrank.eq.0) THEN
 !     call system_clock(i,j,k)
 !     time_end=real(i,kind=8)/real(j,kind=8)
 !     time_cpu = time_end-time_begin
 !     write(*,*)'reduce time=',time_cpu
!     ENDIF
    
!DEC$ENDIF
     ENDIF
        
! ALUG
    IF(myrank.EQ.0) THEN
        
      IF(nlv_glo.EQ.0) THEN
          CALL STIFF_EXACT(i_dir,ipar,nnodeG,nnzG,iaG,jaG,juG,auG)
      ELSE
          CALL STIFF_COARSE(i_dir,ipar,nnodeG,nnzG,iaG,jaG,juG,auG)
      ENDIF
!         
    ENDIF
!
      RETURN
    END
! - - - - - - - - - 
    
 SUBROUTINE stiffness_GC_all(id,nintf,nnode,nnz,ia,au)
          
      USE MD_MG_Global_C, ONLY: i_dir,nlv_glo,nnodeG,nnzG,eG,rG,rG0,imapG,iaG,jaG,juG,auG,auG0,imapGZ, &
                                igather, nsengatA, irevgatA, idispA , imapgatA
      USE MD_MPI, ONLY: myrank
      USE MD_parameter, ONLY: ndom,ipar
      use omp_lib
      
      IMPLICIT NONE
          
!DEC$IF defined (mpi_flag)
      INCLUDE 'mpif.h'
!DEC$ENDIF
!
! input: 
  INTEGER nintf,nnode,nnz,id
  INTEGER ia(nnode+1)
  REAL*8 au(nnz)
! out:
! temp:
  INTEGER i,j,ierr,i1,i2,k,nd,ip

  REAL*8 time_begin,time_end,time_cpu 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
!step-1: transfer residual from local to global
  
   !$omp PARALLEL DO
     DO i=1,nnzG
     auG(i) = 0.d0
     auG0(i) = 0.d0
	 ENDDO
   !$omp end PARALLEL DO
         
    IF(igather.EQ.0) THEN
        
         
         DO i = 1,nintf
            i1 = ia(i)
                i2 = ia(i+1)
                 
            DO j=i1,i2-1
                auG0(imapGZ(j)) = au(j)
                
            ENDDO
        
            
         ENDDO   
         

    ENDIF
    
!          
         IF(ndom.EQ.1) THEN
             
             IF(igather.EQ.0) THEN
             auG = auG0
             ELSE
             !$omp PARALLEL DO
             DO i=1,nnzG			 
             auG(imapgatA(i)) = au(i)
			 ENDDO
             !$omp end PARALLEL DO
             ENDIF
             
         ELSE
             
!DEC$IF defined (mpi_flag)
!    call mpi_barrier(mpi_comm_world,ierr)

! step-2: S&R to all processors

!      call system_clock(i,j,k)
!      time_begin=real(i,kind=8)/real(j,kind=8) 
    IF(igather.EQ.0) THEN
        
       CALL MPI_ALLREDUCE(auG0,auG,nnzG,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
!        CALL MPI_REDUCE(auG0,auG,nnzG,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)

        ELSE
    
              
    CALL MPI_ALLGATHERV(au,nsengatA,mpi_double_precision,auG0,irevgatA,idispA,mpi_double_precision,mpi_comm_world, ierr)
                
    !$omp PARALLEL DO private(i)
            DO j = 1,nnzG 
            i = imapgatA(j) 
            auG(i) = auG0(j)
            ENDDO
    !$omp end PARALLEL DO
                    
         ENDIF

!    call mpi_barrier(mpi_comm_world,ierr)

!      IF(myrank.EQ.0) THEN
!      call system_clock(i,j,k)
!      time_end=real(i,kind=8)/real(j,kind=8)
!      time_cpu = time_end-time_begin
!      write(*,*)'reduce time=',time_cpu

!     ENDIF
    
!DEC$ENDIF
         ENDIF
        
! ALUG
!    IF(myrank.EQ.0) THEN
        
      IF(nlv_glo.EQ.0) THEN
          CALL STIFF_EXACT(i_dir,ipar,nnodeG,nnzG,iaG,jaG,juG,auG)
      ELSE
          IF(id.EQ.1) THEN
          CALL STIFF_COARSE(i_dir,ipar,nnodeG,nnzG,iaG,jaG,juG,auG)
          ELSE
          CALL STIFF_COARSE2(i_dir,ipar,nnodeG,nnzG,iaG,jaG,juG,auG)
          ENDIF
          
      ENDIF
!         
!    ENDIF
!
      RETURN
      END
! - - - - - - - - - 
! = = = = = = = = = = = = = = = = = = = = = = = = = = = !
    SUBROUTINE STIFF_EXACT(i_dir,ipar,nnode,nnz,ia,ja,ju,au)
!
    USE MD_MG_Global_C, ONLY: Ainv, aluG
    use omp_lib
!
    IMPLICIT NONE 
! 
    INTEGER (4) i_dir,ipar
	INTEGER (4) nnode,nnz
	INTEGER (4) ia(*),ja(*),ju(*)
    REAL(8)  au(*)
! 
    INTEGER(4) i,j,i1,i2

       IF(i_dir.NE.0) THEN
           
      !$omp PARALLEL DO private(i1,i2,j)
	  
           DO i=1,nnode
		      Ainv(i,1:nnode) = 0.d0
               i1 = ia(i)
               i2 = ia(i+1)-1
               DO j=i1,i2
                   Ainv(i,ja(j)) = au(j)
               ENDDO
           ENDDO
		   
      !$omp end PARALLEL DO
!           
         IF(nnode.LE.100) THEN
          CALL matrix_inverse_GS_n(nnode,Ainv)
         ELSE
          CALL linalg_invM(nnode, Ainv, 1)
         ENDIF
         
       ENDIF
!         
!
      RETURN
      END
	  
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
! Galerkin Formula: 
    
SUBROUTINE stiff_coarse_global(nnode,nnode1,nfmax,nnz,                                &
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
      integer nnode,nnode1,nfmax
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
      allocate(vi(nnode))
      
    !$omp PARALLEL
    !$omp DO
      DO i=1,nnz1 	
      au1(i) = 0.d0
	  ENDDO
	!$omp END DO
	
!      ni = 0
    !$omp DO
      DO i=1,nnode
      vi(i) = 0.d0 
      ENDDO
	!$omp END DO	  
!      pi = 0.d0
!      nj = 0
!      pj = 0.d0
    !$omp END PARALLEL 
! ---
    !$omp PARALLEL do private(i1,i2,imax,ni,pi,nj,pj,id,k,l,ll,j,j1,j2,jmax,s) firstprivate(vi)
!    !$omp  DO
      
      DO i = 1, nnode1
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
      
!    !$omp end DO
    !$omp end PARALLEL do
! -------------------------------------------------------------!
      DEALLOCATE(vi)
      
      return
      
    END

! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !

! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
      SUBROUTINE STIFF_COARSE(i_dir,ipar,nnodeG,nnzG,iaG,jaG,juG,auG)
! ---
      USE MD_MG_matrix
	  USE MD_MG_coord, only: nnode1,nnode2
	  USE MD_MG_Global_C,ONLY: nlv_glo,inmaxG,inmaxGC,                              &
	                      nnziG,XrestG,jarG,iarG,                                         &
						  nnodeGC,nnzGC,iaGC,jaGC,juGC,auGC,                              &
						  nnziGC,iarGC,jarGC,XrestGC
      use omp_lib
!
      implicit none  
      
      INTEGER(4) i_dir,ipar
      INTEGER(4) nnodeG,nnzG
      INTEGER(4) iaG(nnodeG+1),jaG(nnzG),juG(nnodeG)
      REAL(8) auG(nnzG)
!     
      integer(4) jmax,i,j,k,ilv,nnzr2
      REAL(8) tmp
      real*8 time_begin,time_end,time_cpu
      real*8 time_begin0,time_end0,time_cpu0
	  
!    
      
!      call system_clock(i,j,k)
!      time_begin0=real(i,kind=8)/real(j,kind=8)  
	  
! initial fine level
      nnode1 = nnodeG
      nnz1 = nnzG
      ALLOCATE(ia1(nnode1+1),ja1(nnz1),au1(nnz1))
      
    !$omp PARALLEL
	
    !$omp DO
      DO i=1,nnode1+1
      ia1(i) = iaG(i)
	  ENDDO
    !$omp END DO	  
	  
    !$omp DO
      DO i=1,nnz1
      ja1(i) = jaG(i)
      au1(i) = auG(i)
      ENDDO
    !$omp END DO
	
    !$omp END PARALLEL
      
! 
      DO ilv = 1, nlv_glo
          
! coarse level

      nnode2 = nnodeGC(ilv)
      nnz2 = nnzGC(ilv)
      
      ALLOCATE(ia2(nnode2+1),ja2(nnz2),au2(nnz2))
      
    !$omp PARALLEL
	
    !$omp DO
	  DO i=1,nnode2+1
      ia2(i)=iaGC(i, ilv)
	  ENDDO
    !$omp END DO	  
	  
    !$omp DO
	  DO i=1,nnz2	  
      ja2 (i) = jaGC(i, ilv)
	  ENDDO
    !$omp END DO
	
    !$omp END PARALLEL 
      
! for Xrest:

      IF(ilv.EQ.1) THEN
	  nnzr2 = nnziG
      ALLOCATE(iar2(nnode2+1),jar2(nnzr2),Xrest2(nnzr2))
	  
    !$omp PARALLEL
	
    !$omp DO
	  DO i=1, nnode2+1
      iar2(i)=iarG (i)
	  ENDDO
    !$omp END DO	
	
    !$omp DO
	  DO i=1, nnzr2
      jar2 (i)=jarG (i)
      Xrest2 (i)=XrestG (i)
	  ENDDO
    !$omp END DO
	
    !$omp END PARALLEL
	   
	  ELSE
      nnzr2 = nnziGC(ilv)
     
      ALLOCATE(iar2(nnode2+1),jar2(nnzr2),Xrest2(nnzr2))

    !$omp PARALLEL
	
    !$omp DO
	  DO i=1, nnode2+1
       iar2(i)=iarGC(i, ilv)
	  ENDDO
    !$omp END DO	  

    !$omp DO
	  DO i=1, nnzr2
       jar2 (i)=jarGC (i, ilv)
       Xrest2 (i)=XrestGC (i, ilv)
	  ENDDO
    !$omp END DO

    !$omp END PARALLEL
	   
	   ENDIF
!
       jmax = 4*inmaxGC(ilv)
       jmax = MIN(jmax,nnz1)
       IF(jmax.LE.20) jmax = MIN(nnz1,20)
! 
      CALL stiff_coarse_global(nnode1,nnode2,jmax,                           &
                       nnz1,nnz2,nnzr2,ia1,ja1,ia2,ja2,iar2,jar2,au1,au2,Xrest2)     
      
      DEALLOCATE(iar2,jar2,Xrest2)
       
! adding to array
    !$omp PARALLEL DO

	  DO i=1, nnz2
       auGC(i, ilv) = au2 (i)
	  ENDDO

    !$omp END PARALLEL DO
! update for fine level
       DEALLOCATE(ia1,ja1,au1)
       nnode1 = nnode2
       nnz1 = nnz2
       
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
       
      ENDDO
    
! for the coarsest level:
	 
	 ALLOCATE(ju1(nnode1))
	 ju1(1:nnode1) = juGC(1:nnode1,nlv_glo)
	 
     CALL STIFF_EXACT(i_dir,ipar,nnode1,nnz1,ia1,ja1,ju1,au1)
      
      DEALLOCATE(ia1,ja1,ju1,au1)   
! 
!      call system_clock(i,j,k)
!      time_end0=real(i,kind=8)/real(j,kind=8)
!	  time_cpu0 = time_end0-time_begin0
!      write(*,*)'cpu Gar coarse=',time_cpu0 
    
!     write(*,*)'Gar-post'
      return
	  
    End
    
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
! Galerkin Formula: 
    
SUBROUTINE stiff_coarse_global2(nnode,nnode1,nfmax,nnz,                                &
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
      implicit none
! ---inlet
      integer nnode,nnode1,nfmax
      integer nnz,nnz1,nnzr
      integer ia(*),ja(*)
      integer ia1(*),ja1(*)
      integer iar(*),jar(*)
      real*8 au(*),Xr(*)
!      
! ---outlet 
      real*8 au1(*)
! ---temp
      real*8 vi(nnode)      
      integer i,j,k,l,i1,i2,j1,j2,id,imax,jmax,ll
      real*8 s
      INTEGER(4) ni(nfmax),nj(nfmax)
      REAL(8) pi(nfmax),pj(nfmax)
!      real*8,dimension(:),allocatable::vi
      
! ----------------------------------------------------------!
!      allocate(vi(nnode))
      
      au1(1:nnz1) = 0.d0
      ni = 0
      vi(1:nnode) = 0.d0 
      pi = 0.d0
      nj = 0
      pj = 0.d0
! ---
      DO i = 1, nnode1
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
! -------------------------------------------------------------!
!      DEALLOCATE(vi)
      
      return
      
    END

! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
      SUBROUTINE STIFF_COARSE2(i_dir,ipar,nnodeG,nnzG,iaG,jaG,juG,auG)
! ---
      USE MD_MG_matrix, ONLY: nnz1,nnz2      !, ia1,ja1,ju1,au1, ia2,ja2,au2, iar2, jar2, Xrest2
	  USE MD_MG_coord, only: nnode1,nnode2
	  USE MD_MG_Global_C,ONLY: nlv_glo,inmaxG,inmaxGC,                              &
	                      nnziG,XrestG,jarG,iarG,                                         &
						  nnodeGC,nnzGC,iaGC,jaGC,juGC,auGC,                              &
						  nnziGC,iarGC,jarGC,XrestGC
!
      implicit none  
      
      INTEGER(4) i_dir,ipar
      INTEGER(4) nnodeG,nnzG
      INTEGER(4) iaG(nnodeG+1),jaG(nnzG),juG(nnodeG)
      REAL(8) auG(nnzG)
!     
      integer(4) jmax,i,j,k,ilv,nnzr2
      REAL(8) tmp
      real*8 time_begin,time_end,time_cpu
      real*8 time_begin0,time_end0,time_cpu0
	  
!    temp
      INTEGER(4) ia1(nnodeGC(1)+1), ja1(nnzGC(1))
      REAL(8)  au1(nnzGC(1))
      INTEGER(4) ia2(nnodeGC(1)+1), ja2(nnzGC(1))
      REAL(8)  au2(nnzGC(1))
      INTEGER(4) iar2(nnodeGC(1)+1), jar2(nnziG)
      REAL(8)  Xrest2(nnziG)      
      INTEGER(4) ju1(nnodeGC(nlv_glo))
      
!      call system_clock(i,j,k)
!      time_begin0=real(i,kind=8)/real(j,kind=8)  
	  
! initial fine level
      nnode1 = nnodeG
      nnz1 = nnzG
!      ALLOCATE(ia1(nnode1+1),ja1(nnz1),au1(nnz1))
!      ia1 = iaG
!      ja1 = jaG
!      au1 = auG
! 
      DO ilv = 1, nlv_glo
          
! coarse level

      nnode2 = nnodeGC(ilv)
      nnz2 = nnzGC(ilv)
      
!      ALLOCATE(ia2(nnode2+1),ja2(nnz2),au2(nnz2))
      
      ia2(1:nnode2+1)=iaGC(1:nnode2+1, ilv)
      ja2 (1:nnz2) = jaGC(1:nnz2, ilv)
      
! for Xrest:

      IF(ilv.EQ.1) THEN
	  nnzr2 = nnziG
!      ALLOCATE(iar2(nnode2+1),jar2(nnzr2),Xrest2(nnzr2))
	  
      iar2(1:nnode2+1)=iarG (1:nnode2+1)

      jar2 (1:nnzr2)=jarG (1:nnzr2)

      Xrest2 (1:nnzr2)=XrestG (1:nnzr2)
	   
	  ELSE
      nnzr2 = nnziGC(ilv)
     
!      ALLOCATE(iar2(nnode2+1),jar2(nnzr2),Xrest2(nnzr2))

       iar2(1:nnode2+1)=iarGC(1:nnode2+1, ilv)

       jar2 (1:nnzr2)=jarGC (1:nnzr2, ilv)

       Xrest2 (1:nnzr2)=XrestGC (1:nnzr2, ilv)
	   
	   ENDIF
!
       jmax = 4*inmaxGC(ilv)
       jmax = MIN(jmax,nnz1)
       IF(jmax.LE.20) jmax = MIN(nnz1,20)
! 
       IF(ilv.EQ.1) THEN
         CALL stiff_coarse_global2(nnode1,nnode2,jmax,                           &
                       nnz1,nnz2,nnzr2,iaG,jaG,ia2,ja2,iar2,jar2,auG,au2,Xrest2) 
       ELSE
         CALL stiff_coarse_global2(nnode1,nnode2,jmax,                           &
                       nnz1,nnz2,nnzr2,ia1,ja1,ia2,ja2,iar2,jar2,au1,au2,Xrest2) 
       ENDIF
           
      
!      DEALLOCATE(iar2,jar2,Xrest2)
       
! adding to array
       auGC(1:nnz2, ilv) = au2 (1:nnz2)
! update for fine level
!       DEALLOCATE(ia1,ja1,au1)
       nnode1 = nnode2
       nnz1 = nnz2
       
!       ALLOCATE(ia1(nnode1+1),ja1(nnz1),au1(nnz1))
       ia1(1:nnode1+1) = ia2(1:nnode1+1)
       ja1(1:nnz1) = ja2(1:nnz1)
       au1(1:nnz1) = au2(1:nnz1)
 !      DEALLOCATE(ia2,ja2,au2)
       
      ENDDO
    
! for the coarsest level:
	 
      IF(nnode1.NE.nnodeGC(nlv_glo)) THEN
          WRITE(999,*)'error in Stiffness_GC2,',nnode1,nnodeGC(nlv_glo)
      ENDIF
!
!	 ALLOCATE(ju1(nnode1))
	 ju1(1:nnode1) = juGC(1:nnode1,nlv_glo)
	 
     CALL STIFF_EXACT(i_dir,ipar,nnode1,nnz1,ia1,ja1,ju1,au1)
      
!      DEALLOCATE(ia1,ja1,ju1,au1)   
! 
!      call system_clock(i,j,k)
!      time_end0=real(i,kind=8)/real(j,kind=8)
!	  time_cpu0 = time_end0-time_begin0
!      write(*,*)'cpu Gar coarse=',time_cpu0 
    
!     write(*,*)'Gar-post'
      return
	  
    End