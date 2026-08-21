! - - - - - - - - - 
    
 SUBROUTINE stiffness_GC_all(nintf,nnode,nnz,ia,au)
          
      USE MD_MG_Global_C, ONLY: i_dir,nlv_glo,nnodeG,nnzG,eG,rG,rG0,imapG,iaG,jaG,juG,auG,auG0,imapGZ, &
                                igather, nsengatA, irevgatA, idispA , imapgatA
      USE MD_MPI, ONLY: myrank
      USE MD_parameter, ONLY: ndom
      use omp_lib
      
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

! step-2: S&R to all processors

    IF(igather.EQ.0) THEN
        
       CALL MPI_ALLREDUCE(auG0,auG,nnzG,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)

        ELSE
    
              
    CALL MPI_ALLGATHERV(au,nsengatA,mpi_double_precision,auG0,irevgatA,idispA,mpi_double_precision,mpi_comm_world, ierr)
                
    !$omp PARALLEL DO private(i)
            DO j = 1,nnzG 
            i = imapgatA(j) 
            auG(i) = auG0(j)
            ENDDO
    !$omp end PARALLEL DO
                    
         ENDIF



    
!DEC$ENDIF
         ENDIF
        
! ALUG
        
      IF(nlv_glo.EQ.0) THEN
          CALL STIFF_EXACT(i_dir,nnodeG,nnzG,iaG,jaG,juG,auG)
      ELSE
          CALL STIFF_COARSE2(i_dir,nnodeG,nnzG,iaG,jaG,juG,auG)
      ENDIF
!         
!
      RETURN
      END
! - - - - - - - - - 
! = = = = = = = = = = = = = = = = = = = = = = = = = = = !
    SUBROUTINE STIFF_EXACT(i_dir,nnode,nnz,ia,ja,ju,au)
!
    USE MD_MG_Global_C, ONLY: Ainv, aluG
    use omp_lib
!
    IMPLICIT NONE 
! 
    INTEGER (4) i_dir
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
      RETURN
      END
	  
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
! Galerkin Formula: 
    

! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !

! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
    
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
! Galerkin Formula: 
    
SUBROUTINE stiff_coarse_global2(nnode,nnode1,nfmax,nnz,                                &
                              nnz1,nnzr,ia,ja,ia1,ja1,iar,jar,au,au1,Xr)
    
! ************************************************************************************!
! this subroutine calculates the coarse-stiff matrix by Galerkin formular             !
!   or: Ac(I,J) = R(I,k)*Af(k,l)*R(J,l)                                               !
!   for each I, J on coarse grid, we only find k,l on fine-grid such that:            !
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
      
! ----------------------------------------------------------!
      
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
      
      return
      
    END

! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
      SUBROUTINE STIFF_COARSE2(i_dir,nnodeG,nnzG,iaG,jaG,juG,auG)
! ---
      USE MD_MG_matrix, ONLY: nnz1,nnz2      !, ia1,ja1,ju1,au1, ia2,ja2,au2, iar2, jar2, Xrest2
	  USE MD_MG_coord, only: nnode1,nnode2
	  USE MD_MG_Global_C,ONLY: nlv_glo,inmaxG,inmaxGC,                              &
	                      nnziG,XrestG,jarG,iarG,                                         &
						  nnodeGC,nnzGC,iaGC,jaGC,juGC,auGC,                              &
						  nnziGC,iarGC,jarGC,XrestGC
!
      implicit none  
      
      INTEGER(4) i_dir
      INTEGER(4) nnodeG,nnzG
      INTEGER(4) iaG(nnodeG+1),jaG(nnzG),juG(nnodeG)
      REAL(8) auG(nnzG)
!     
      integer(4) jmax,i,j,k,ilv,nnzr2
      REAL(8) tmp
	  
!    temp
      INTEGER(4) ia1(nnodeGC(1)+1), ja1(nnzGC(1))
      REAL(8)  au1(nnzGC(1))
      INTEGER(4) ia2(nnodeGC(1)+1), ja2(nnzGC(1))
      REAL(8)  au2(nnzGC(1))
      INTEGER(4) iar2(nnodeGC(1)+1), jar2(nnziG)
      REAL(8)  Xrest2(nnziG)      
      INTEGER(4) ju1(nnodeGC(nlv_glo))
      
	  
! initial fine level
      nnode1 = nnodeG
      nnz1 = nnzG
! 
      DO ilv = 1, nlv_glo
          
! coarse level

      nnode2 = nnodeGC(ilv)
      nnz2 = nnzGC(ilv)
      
      
      ia2(1:nnode2+1)=iaGC(1:nnode2+1, ilv)
      ja2 (1:nnz2) = jaGC(1:nnz2, ilv)
      
! for Xrest:

      IF(ilv.EQ.1) THEN
	  nnzr2 = nnziG
	  
      iar2(1:nnode2+1)=iarG (1:nnode2+1)

      jar2 (1:nnzr2)=jarG (1:nnzr2)

      Xrest2 (1:nnzr2)=XrestG (1:nnzr2)
	   
	  ELSE
      nnzr2 = nnziGC(ilv)
     

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
           
      
       
! adding to array
       auGC(1:nnz2, ilv) = au2 (1:nnz2)
! update for fine level
       nnode1 = nnode2
       nnz1 = nnz2
       
       ia1(1:nnode1+1) = ia2(1:nnode1+1)
       ja1(1:nnz1) = ja2(1:nnz1)
       au1(1:nnz1) = au2(1:nnz1)
       
      ENDDO
    
! for the coarsest level:
	 
      IF(nnode1.NE.nnodeGC(nlv_glo)) THEN
          WRITE(*,*)'PMG error in Stiffness_GC2,',nnode1,nnodeGC(nlv_glo)
      ENDIF
!
	 ju1(1:nnode1) = juGC(1:nnode1,nlv_glo)
	 
     CALL STIFF_EXACT(i_dir,nnode1,nnz1,ia1,ja1,ju1,au1)
      
! 
    
      return
	  
    End