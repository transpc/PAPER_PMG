      SUBROUTINE e_neighbor_e_nnz(ncn,nnode,nelem,nnd,inode,nnz,ia,ja,nnz_neigh)

!
      USE MD_connectivity, ONLY:  live
      USE MD_MG_index, ONLY: report_text
      
	  IMPLICIT NONE
	  
! this sub. conducts neighbor element of each element (shared face)   !
! used for metis sub.    											  !
! modify by Dr.Ha (4/2023)	
! ncn: Number Common Codes of two neighbor element. 
! only need for FEM
! 
! ---inlet 
      INTEGER(4) nnd,ncn
      INTEGER(4) nnode,nelem,nnz
      INTEGER(4) inode(nnd,nelem)
      INTEGER(4) ia(nnode+1),ja(nnz)
! out:
      INTEGER(4) nnz_neigh

! --- temp
      INTEGER(4):: alstatus
      INTEGER(4) ID(nnd),ie,nn,ne,I,J,JD,index, J1, J2

! ----------------------------------------!       

! ---    
         Allocate(live(nnode),stat=alstatus)
         
        IF (alstatus/=0) THEN
         report_text = 'not enough memory,serial-e_neighbor_e,live'
         CALL STOP_MPI(report_text)
        ENDIF
        
         live=0  
!
         nn = 0
         
      DO 100 ne=1,nelem

         ID(1:nnd)=inode(1:nnd,ne)
! ---
         live(ID(1:nnd))=1
		 
! loop for nnd-1 vertices of element:

        DO I = 1,nnd-1
		 
		  JD = ID(I)
		 
          J1 = ia(JD)
          J2 = ia(JD+1)-1
		  DO J = J1,J2
		  
             ie=ja(J)      ! neighbor element of JD
			 
             index=SUM(live(inode(1:nnd,ie)))
			 
             IF(index.EQ.ncn)THEN

             nn=nn+1
             
             ENDIF
			 
		  ENDDO
		  
          live(JD)=-100
		  
        ENDDO
		  
          live(ID(1:nnd))=0
! ---
100 ENDDO
    
      nnz_neigh = nn

!
       DEALLOCATE(live)
!
        RETURN 
    END
    
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

    
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
    
    SUBROUTINE e_neighbor_e(ncn,nnode,nelem,nnd,inode,nnz,ia,ja,nnz_neigh,num_neigh,neigh)

!
      USE MD_connectivity, ONLY:  live
      USE MD_MG_index, ONLY: report_text
      
	  IMPLICIT NONE
	  
! this sub. conducts neighbor element of each element (shared face)   !
! used for metis sub.    											  !
! modify by Dr.Ha (4/2023)	
! ncn: Number Common Codes of two neighbor element. 
! only need for FEM
! 
! ---inlet 
      INTEGER(4) nnd,ncn
      INTEGER(4) nnode,nelem,nnz,nnz_neigh
      INTEGER(4) inode(nnd,nelem)
      INTEGER(4) ia(nnode+1),ja(nnz)
! out:
      INTEGER(4) neigh(nnz_neigh),num_neigh(nelem+1)

! --- temp
      INTEGER(4):: alstatus
      INTEGER(4) ID(nnd),ie,nn,ne,I,J,JD,index, J1, J2

! ----------------------------------------!       

! ---    
         Allocate(live(nnode),stat=alstatus)
         
        IF (alstatus/=0) THEN
         report_text = 'not enough memory,serial-e_neighbor_e,live'
         CALL STOP_MPI(report_text)
        ENDIF
        
         live=0  
!
         num_neigh(1) = 1
         
      DO 100 ne=1,nelem

         ID(1:nnd)=inode(1:nnd,ne)
! ---
         live(ID(1:nnd))=1
         
         nn=num_neigh(ne)
		 
! loop for nnd-1 vertices of element:

        DO I = 1,nnd-1
		 
		  JD = ID(I)
		 
          J1 = ia(JD)
          J2 = ia(JD+1)-1
		  DO J = J1,J2
		  
             ie=ja(J)      ! neighbor element of JD
			 
             index=SUM(live(inode(1:nnd,ie)))
			 
             IF(index.EQ.ncn)THEN

             neigh(nn) = ie
             nn=nn+1
             
             ENDIF
			 
		  ENDDO
		  
          live(JD)=-100
		  
        ENDDO
		 
		  num_neigh(ne+1) = nn
		  
          live(ID(1:nnd))=0
! ---
100   ENDDO

! check:
         IF(num_neigh(nelem+1).NE.(nnz_neigh+1)) THEN
           report_text = 'num_neigh(nelem+1).NE.(nnz_neigh+1)'
           CALL STOP_MPI(report_text)
		 ENDIF
!
       DEALLOCATE(live)
!
        RETURN 
    END
    
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

      subroutine e_neighbor_node(nelem,nnode,nnd,inode,mxnbne,iwk,iwork,ierr)
! --- P1P1 element neighbour of vertex
      IMPLICIT NONE
! ---inlet 
      integer(4) nnode,nelem,mxnbne,nnd
      integer(4) inode(nnd,nelem)
! --- out 
      INTEGER(4) ierr
      integer(4) iwk(nnode),iwork(mxnbne,nnode)
!    temp
      integer(4) ne,j,n,k
!..
        ierr = 0
! ---
        do ne=1,nelem 
           do j=1,nnd 
              n = inode(j,ne)
              k = iwk(n) + 1
              IF(k.GT.mxnbne) THEN
                  ierr = 1
                  RETURN
              ENDIF
               
              iwk(n) = k
              iwork(k,n) = ne
           enddo
        enddo

        return 
    end 
    
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    subroutine node_neighbor_node_FEM_nnz(nvpe,nnode,nelem,inode,nnz_en, ia_en, ja_en,nnz_neigh,nf)
	  
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * !
! ---
      implicit none
! ---inlet 
      integer(4) nnode,nelem,nvpe, nnz_en
      integer(4) inode(nvpe,nelem)
      INTEGER(4):: ia_en(nnode+1), ja_en(nnz_en)
! --- out 
      integer(4) nnz_neigh,nf
!    temporary
      INTEGER(4):: alstatus
      integer(4) i,j,k,ne,n,nnd,id,nd, j1, j2
      integer(4),dimension(:),allocatable:: imark
!    		
      allocate(imark(nnode),stat=alstatus)
        imark = 0
      
        IF (alstatus/=0) THEN
          WRITE(*,*)'not enough memory,serial-node-neighbor_node'
          STOP
        ENDIF
!	   
        nnd = 0
        nf = 0
        DO i = 1,nnode
         
         imark(i)=i
         j1 = ia_en(i)
         j2 = ia_en(i+1)-1
         
         nd = 0
         do j = j1,j2
           ne = ja_en(j)
           do n = 1,nvpe
		    id = inode(n,ne)
			
            if(imark(id).ne.i) then
            nnd = nnd + 1
            nd = nd+1
            
            imark(id) = i
            end if
			
          end do
         end do
         
		 nf = MAX(nf, nd)
         
        ENDDO

        nnz_neigh = nnd
! ---
   	  Deallocate(imark)
	  
      return
    END
    
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
    
    subroutine node_neighbor_node_FEM(nvpe,nnode,nelem,inode,nnz_en, ia_en, ja_en,nnz_neigh, ia_neigh,ja_neigh)
	  
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * !
! ---
      implicit none
! ---inlet 
      integer(4) nnode,nelem,nvpe, nnz_en, nnz_neigh
      integer(4) inode(nvpe,nelem)
      INTEGER(4):: ia_en(nnode+1), ja_en(nnz_en)
! --- out 
      integer(4) ia_neigh(nnode+1),ja_neigh(nnz_neigh)
!    temporary
      INTEGER(4):: alstatus
      integer(4) i,j,k,ne,n,nnd,id,nd, j1, j2
      integer(4),dimension(:),allocatable:: imark
!    		
      allocate(imark(nnode),stat=alstatus)
        imark = 0
      
        IF (alstatus/=0) THEN
          WRITE(*,*)'not enough memory,serial-node-neighbor_node'
          STOP
        ENDIF
!	   
        ia_neigh(1) = 1
!    
        DO i = 1,nnode
		 nnd = ia_neigh(i)-1
         
         imark(i)=i
         j1 = ia_en(i)
         j2 = ia_en(i+1)-1
         
         do j = j1,j2
           ne = ja_en(j)
           do n = 1,nvpe
		    id = inode(n,ne)
			
            if(imark(id).ne.i) then
            nnd = nnd + 1
            ja_neigh(nnd) = id
            imark(id) = i
            end if
			
          end do
         end do
         
         
		 
         ia_neigh(i+1) = nnd+1
		 
        ENDDO

! ---
   	  Deallocate(imark)
	  
      return
    END
    
! - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - - - - 
      SUBROUTINE neighbor_fine_graph(nnodef, nnodec, nnz_neighc, ia_neighc, ja_neighc,  &
                                        icoarse, nnz,jmax, ia, ja)
!
         IMPLICIT NONE
!
!.....input
!
         INTEGER(4) nnodec, nnodef, jmax, nnz, nnz_neighc
         INTEGER(4) ia_neighc(nnodef+1), ja_neighc(nnz_neighc), icoarse(nnodef)

!
!.....output
!
         INTEGER(4) ia(nnodef+1), ja(nnz)
!
!.....temporary
!
         INTEGER(4) i, j, k, l, nnd, ie, je,nd, nn, i1, i2
         INTEGER(4) ni(jmax)
!
         ni = 0
!
         ia(1) = 1
         
         DO ie = 1, nnodef

!    
            nd = ia(ie)
            
            k = icoarse(ie)
            IF (k .NE. 0) THEN
            ja(nd) = k
            ia(ie+1) = nd+1
            CYCLE
            END IF

!

            i1 = ia_neighc(ie)
            i2 = ia_neighc(ie+1)-1
            nn = i2-i1+1
            
            nnd = 0
            
            DO i = 1, nn                   !nnei(ie)
               je = ja_neighc(i1+i-1)
               k = icoarse(je)
               IF (k .EQ. 0) CYCLE
               nnd = nnd+1
               ni(nnd) = k
            END DO
!

! ordering
! it will be re-ordered later after reducing no. ni.
            
!
            DO i = 1, nnd
            
            ja(nd) = ni(i)
            nd = nd + 1
            ENDDO
            
            ia(ie+1) = nd

         END DO

! -------------------------------------!
         RETURN
    END
! - - - - - - - - - - - - - - - - - - - !
    
! - - - - - - - - - - - - - - - - - - - - - - - !
      !
      SUBROUTINE neighbor_fine_graph2(nnodef, nnodec, nnz_neighc, ia_neighc, ja_neighc,  &
                                        icoarse, nnz,jmax, ia, ja)
!
         IMPLICIT NONE
!
!.....input
!
         INTEGER(4) nnodec, nnodef, jmax, nnz, nnz_neighc
         INTEGER(4) ia_neighc(nnodef+1), ja_neighc(nnz_neighc), icoarse(nnodef)

!
!.....output
!
         INTEGER(4) ia(nnodef+1), ja(nnz)
!
!.....temporary
!
         INTEGER(4) i, j, k, l, nnd, ie, je,nd, nn, i1, i2, j1, j2, jj, nn2
         INTEGER(4) ni(jmax)
         INTEGER,DIMENSION(:),ALLOCATABLE :: imark
!
         allocate(imark(nnodec))
         imark = 0
         
         ni = 0
!
         ia(1) = 1
         
         DO ie = 1, nnodef

!    
            nd = ia(ie)
            
            k = icoarse(ie)
            IF (k .NE. 0) THEN
            ja(nd) = k
            ia(ie+1) = nd+1
            CYCLE
            END IF

!
            i1 = ia_neighc(ie)
            i2 = ia_neighc(ie+1)-1
            nn = i2-i1+1
            
            nnd = 0
            
            DO i = 1, nn                   !nnei(ie)
               je = ja_neighc(i1+i-1)
               k = icoarse(je)
              IF (k .NE. 0) THEN
                   
                 IF(imark(k).NE.ie) THEN
                 nnd = nnd+1
                 ni(nnd) = k
                 imark(k) = ie
                 ENDIF
               
               ELSE
                   
                 j1 = ia_neighc(je)
                 j2 = ia_neighc(je+1)-1
                 nn2 = j2-j1+1
                 
                 DO j=1,nn2
                  jj = ja_neighc(j1+j-1)
                  k = icoarse(jj)
                    IF(k.NE.0) THEN
                      IF(imark(k).NE.ie) THEN
                      nnd = nnd+1
                      ni(nnd) = k
                      imark(k) = ie
                      ENDIF
                    ENDIF
                 ENDDO
                 
               ENDIF
!  
               
            END DO
!

! ordering
! ordering
! it will be re-ordered later after reducing no. ni.
            
!
            DO i = 1, nnd
            
            ja(nd) = ni(i)
            nd = nd + 1
            ENDDO
            
            ia(ie+1) = nd

         END DO

      DEALLOCATE(imark)	
! -------------------------------------!
         RETURN
    END
    
! - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - - - - - - - - - - - - - - - - - - !
          !
      SUBROUTINE neighbor_fine_graph_nnz(nnodef, nnodec, nnz_neighc, ia_neighc, ja_neighc,  &
                 icoarse, nnz, jmax)
!
         USE MD_MPI, ONLY: myrank
!/
         IMPLICIT NONE
!
!.....input
!
         INTEGER(4) nnodec, nnodef, nnz_neighc
         INTEGER(4) ia_neighc(nnodef+1), ja_neighc(nnz_neighc), icoarse(nnodef)

!
!.....output
!
         INTEGER(4)  nnz, jmax
!
!.....temporary
!
         INTEGER(4) i, j, k, l, nnd, ie, je,nd, nn, i1, i2
!
         nnz = 0
         jmax = 1
         
         DO ie = 1, nnodef
!       
            k = icoarse(ie)
            IF (k .NE. 0) THEN
            nnz = nnz +1
            CYCLE
            END IF

!

            i1 = ia_neighc(ie)
            i2 = ia_neighc(ie+1)-1
            nn = i2-i1+1
            
            nnd = 0
            
            DO i = 1, nn                   !nnei(ie)
               je = ja_neighc(i1+i-1)
               k = icoarse(je)
               IF (k .EQ. 0) CYCLE
               nnd = nnd+1
            END DO
!  
            
            nnz = nnz + nnd
            
            jmax = MAX(jmax, nnd)

         END DO
         
 IF(myrank == 0) THEN
    WRITE(999,*)'nnz and jmax for interpolation', nnz, jmax
 ENDIF
 

! -------------------------------------!
         RETURN
    END
! - - - - - - - - - - - - - - - - - - - !
    
! - - - - - - - - - - - - - - - - - - - - - - - !
      !
      SUBROUTINE neighbor_fine_graph2_nnz(nnodef, nnodec, nnz_neighc, ia_neighc, ja_neighc,  &
                                        icoarse, nnz, jmax)
!
         IMPLICIT NONE
!
!.....input
!
         INTEGER(4) nnodec, nnodef, nnz_neighc
         INTEGER(4) ia_neighc(nnodef+1), ja_neighc(nnz_neighc), icoarse(nnodef)

!
!.....output
!
         INTEGER(4) nnz, jmax
!
!.....temporary
!
         INTEGER(4) i, j, k, l, nnd, ie, je,nd, nn, i1, i2, j1, j2, jj, nn2
         INTEGER,DIMENSION(:),ALLOCATABLE :: imark
!
         allocate(imark(nnodec))
         imark = 0
!
         
         nnz = 0
         jmax = 1
         
         DO ie = 1, nnodef

!    
            k = icoarse(ie)
            IF (k .NE. 0) THEN
            nnz = nnz + 1
            CYCLE
            END IF

!
            i1 = ia_neighc(ie)
            i2 = ia_neighc(ie+1)-1
            nn = i2-i1+1
            
            nnd = 0
            
            DO i = 1, nn                   !nnei(ie)
               je = ja_neighc(i1+i-1)
               k = icoarse(je)
              IF (k .NE. 0) THEN
                   
                 IF(imark(k).NE.ie) THEN
                 nnd = nnd+1
                 imark(k) = ie
                 ENDIF
               
               ELSE
                   
                 j1 = ia_neighc(je)
                 j2 = ia_neighc(je+1)-1
                 nn2 = j2-j1+1
                 
                 DO j=1,nn2
                  jj = ja_neighc(j1+j-1)
                  k = icoarse(jj)
                    IF(k.NE.0) THEN
                      IF(imark(k).NE.ie) THEN
                      nnd = nnd+1
                      imark(k) = ie
                      ENDIF
                    ENDIF
                 ENDDO
                 
               ENDIF
!  
               
            END DO
!

            nnz = nnz + nnd
            
            jmax = MAX(jmax, nnd)
            

         END DO

      DEALLOCATE(imark)	
      
    WRITE(999,*)'nnz and jmax for interpolation', nnz, jmax      
! -------------------------------------!
         RETURN
    END


! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
    
      SUBROUTINE reduce_neibor(ndim,jmax,ip_nmax,teta_p,nnodf,nnodc,coordf,coordc,nnz,ia,ja)
!
      IMPLICIT NONE
!  - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
! reduce the number of coarse cell neighbour of cell.                        !
! using distance of fine-central-cell and coarse-cell						 !
! max_nei: maximum number of cell neighbor. 								 !
! nnode:fine; nnode1: coarse 												 !
! designed: shaha - June -2018												 !
! updated: april 2023:                                                       !
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - !
! input
      INTEGER (4) jmax,ip_nmax,ndim, nnz, nnodf, nnodc
      INTEGER (4) ia(nnodf+1), ja(nnz)
      REAL(8) coordf(ndim,nnodf),coordc(ndim,nnodc)
      REAL(8) teta_p
! output (ia,ja)

! free 
      integer i,j,k,l,i1,i2,k1,k2,l1,l2,j1,j2,ie,ne
      integer nnd,id,id1,id2,imax,kmax,nnz_tmp,nnd2
      INTEGER(4) ni(jmax)
      real*8 xtmp, dmin
      integer,dimension(:),allocatable::ia_tmp
      real*8 dx(jmax),xc(3),dx2(jmax)
! --------------------------------------------------!    
      
      allocate(ia_tmp(nnodf+1))
      
      ia_tmp = 0
      dx=0.0d0
      ni = 0

! ---	  
      ia_tmp(1) = 1
      nnz_tmp = 0
!
      
       DO ie = 1,nnodf
	   
           i1 = ia(ie)
           i2 = ia(ie+1)-1
           kmax = i2-i1+1
! 
           IF(kmax.EQ.1) THEN
               nnz_tmp = nnz_tmp + 1
               ja(nnz_tmp) = ja(i1)
               ia_tmp(ie+1) = ia_tmp(ie)+1
               CYCLE
           ENDIF
! kmax>1:
! 
           IF((teta_p.LT.0.1).AND.(kmax.LE.ip_nmax)) THEN                 ! not changed
              DO i = 1,kmax
               nnz_tmp = nnz_tmp + 1
               ja(nnz_tmp) = ja(i1+i-1)
              ENDDO
              
               ia_tmp(ie+1) = ia_tmp(ie)+kmax
               CYCLE              
           ENDIF
       
! need to changed:
	
	      xc(1:ndim) = coordf(1:ndim,ie)
!		  
          dmin = 1.d10
          DO i = 1,kmax
            
            j=ja(i1+i-1)
            ni(i) = j
            
            IF(ndim.EQ.2) THEN
		    dx(i) = (coordc(1,j) - xc(1))**2.d0+(coordc(2,j)-xc(2))**2.d0
            ELSE
		    dx(i) = (coordc(1,j) - xc(1))**2.d0+(coordc(2,j)-xc(2))**2.d0+(coordc(3,j)-xc(3))**2.d0
            ENDIF
                
			dx(i)=DSQRT(dx(i))
            dmin = min(dmin, dx(i))
          ENDDO
		  
! selecting using teta_p
          
          IF(teta_p.GE.0.1) THEN
            dmin = 1.d0/teta_p*dmin
            nnd= 0
            
            DO i = 1, kmax
               xtmp = dx(i)
               IF (xtmp .LE. dmin) THEN
                  nnd = nnd+1
                  ni(nnd) = ni(i)                          ! using only 1 array: ni
                  dx2(nnd) = dx(i)
                  
               END IF
!
            END DO
            
          ELSE
              nnd = kmax 
              dx2(1:nnd) = dx(1:nnd)
          ENDIF
          
! change by ip_nmax
          
        IF(nnd.GT.ip_nmax) THEN
            
          nnd2 = 0
          DO i = 1,nnd
		    xtmp = dx2(i)
		    id = 0
            DO j = 1,nnd
            IF(j.eq.i) CYCLE
            IF(xtmp.GT.dx2(j)) THEN
			id = id+1
            ENDIF
            ENDDO
! 
            IF(id.LE.(ip_nmax-1)) THEN
            nnd2 = nnd2+1
            ni(nnd2) = ni(i)
            ENDIF
!
          ENDDO
          
        ELSE
            nnd2 = nnd
        ENDIF

! update array:
        DO i = 1,nnd2
           nnz_tmp = nnz_tmp + 1
           ja(nnz_tmp) = ni(i)      !ja(i1+i-1)
        ENDDO
              
        ia_tmp(ie+1) = ia_tmp(ie)+nnd2        
        
! 
        dx(1:kmax) = 0.d0
        ni(1:kmax) = 0
        dx2(1:nnd) = 0.d0
        
! 
        
       END DO
       
! Checking 
       
       IF((nnz_tmp+1).NE.ia_tmp(nnodf+1)) THEN
           WRITE(*,*)'error in reduce neighbor',nnz+1, ia_tmp(nnodf+1)
           STOP
       ENDIF
! 
       ia(1:nnodf+1) = ia_tmp(1:nnodf+1)
       
       
      DEALLOCATE(ia_tmp)
! -------------------------------------------------!      
      RETURN
      
      END 

! = = = = = = = = = = = = = = = = = = = = = = = = = !
    