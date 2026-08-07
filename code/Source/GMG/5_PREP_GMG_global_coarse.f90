      SUBROUTINE PREP_GMG_global_coarse

         USE MD_parameter, ONLY: ndim, teta,teta_p, alpha, isemi,ip_inter
         USE MD_geometry, ONLY: imap
         USE MD_MG_coord
         USE MD_MG_index, ONLY:  ip_nmax, ioplv,ip_lev
         USE MD_MG_matrix
		 USE MD_MG_Global_C
         USE MD_connectivity, ONLY: ia_neighc, ja_neighc, nnz_neighc,            &
                                     nnz_tmp, ia_tmp, ja_tmp
! ---
         IMPLICIT NONE

! ---temp
         INTEGER(4) i, k, nv, nd, j, k1, k2, l, i1, i2, nnd
         INTEGER(4) ilv, nnzt
         INTEGER(4) jmax, jmax1

!         INTEGER(4), DIMENSION(:), ALLOCATABLE :: nneit
!         INTEGER(4), DIMENSION(:, :), ALLOCATABLE :: ineit
!         INTEGER(4), DIMENSION(:), ALLOCATABLE :: iwk, iat, jat
!         INTEGER(4), DIMENSION(:, :), ALLOCATABLE :: iwork
!         REAL(8), DIMENSION(:), ALLOCATABLE :: aut

!
         nnode1 = nnodeG
         nnz1 = nnzG
		 
         ALLOCATE (coord1(ndim, nnode1))
         coord1 = coordG
         ALLOCATE (ia1(nnode1+1), ja1(nnz1))
         ia1 = iaG
         ja1 = jaG
		 
         DEALLOCATE(coordG)
!       For level aG - - - - - - - - - - - - - - - - - - 
         IF(nlv_glo.NE.0) THEN
          ALLOCATE(inmaxGC(nlv_glo))
         ENDIF
         

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
!       

         DO ilv = 1, nlv_glomax

! 1: neighbor nodes of each nodes  - - - - - - - - - - !
             
            nnz_neighc = nnz1 - nnode1                         ! only for neighbor nodes
            
            ALLOCATE (ia_neighc(nnode1+1), ja_neighc(nnz_neighc))

            nmax1 = 1
            ia_neighc(1) = 1
            nnzt = 0
            DO i = 1, nnode1
               i1 = ia1(i)
               i2 = ia1(i+1)-1

               nnd = 0
               DO j = i1, i2
                  k = ja1(j)
                  IF (k .EQ. i) CYCLE
                  nnd = nnd+1
                  nnzt = nnzt+1
                  ja_neighc(nnzt) = k
               END DO
               ia_neighc(i+1) = ia_neighc(i)+ nnd
               nmax1 = MAX(nmax1, nnd)

            END DO
            
            IF(ilv.EQ.1) THEN
            inmaxG = nmax1
            ENDIF
            
            inmaxGC(ilv) = nmax1
            

! 2: coarsening step

            nnode2 = INT(0.6*nnode1)
            IF (nnode1 .LE. 2000) nnode2 = nnode1

            ALLOCATE (imap(nnode2), icoarse(nnode1))
            
            imap = 0
            icoarse = 0

               IF (isemi .eq. 0) THEN
                  CALL coarsening_semi(ndim, nnode1, nmax1, nnz_neighc, ia_neighc, ja_neighc, nnode2, imap, icoarse, teta, coord1)
!                  CALL coarsening_semi(ndim, nnode1, nmax1, nnei, inei, nnode2, imap, icoarse, teta, coord1)
               ELSE
 !                 CALL coarsening_semi_amg(ndim, nnode1, nmax1, nnei, inei, nnode2, imap, icoarse, teta, coord1)
               END IF
!
            ALLOCATE (coord2(ndim, nnode2))
!
            DO i = 1, nnode2
               j = imap(i)

               coord2(1:ndim, i) = coord1(1:ndim, j)

            END DO
! 3: Iterpolation procedure
!            jmax = INT(1.5*nmax1)
!            IF(ilv.GE.3) jmax = 3*nmax1
!            IF(ilv.GE.5) jmax = 5*nmax1
!    3.1: iwk-neighbor nodes of each fine-cell

!            ALLOCATE (iwk(nnode1),iwork(jmax,nnode1))
              
              IF(ip_lev.EQ.1) THEN 
               CALL neighbor_fine_graph_nnz(nnode1,nnode2,nnz_neighc, ia_neighc, ja_neighc,icoarse,nnz_tmp, jmax)
              ELSE
               CALL neighbor_fine_graph2_nnz(nnode1,nnode2,nnz_neighc, ia_neighc, ja_neighc,icoarse,nnz_tmp, jmax)
              ENDIF

!
! allocate:
            ALLOCATE(ia_tmp(nnode1+1), ja_tmp(nnz_tmp))
!
              IF(ip_lev.EQ.1) THEN 
               CALL neighbor_fine_graph(nnode1,nnode2,nnz_neighc, ia_neighc, ja_neighc,icoarse,nnz_tmp, jmax, ia_tmp, ja_tmp)
              ELSE
               CALL neighbor_fine_graph2(nnode1,nnode2,nnz_neighc, ia_neighc, ja_neighc,icoarse,nnz_tmp, jmax, ia_tmp, ja_tmp)
              ENDIF
!
            DEALLOCATE(ia_neighc, ja_neighc)
!    3.2: reduce to n-max nodes by order of distance

!    3.2: reduce iwk (using teta)

! notes this one: using reduce by teta or nmax?
        i = 0
        IF(teta_p.GT.0.1) i = i + 1 
        IF((ip_nmax.NE.0).AND.(jmax.GT.ip_nmax)) i = i + 1     
!
        IF(i.NE.0) THEN
        CALL reduce_neibor(ndim,jmax,ip_nmax,teta_p,nnode1,nnode2,coord1,coord2,nnz_tmp,ia_tmp,ja_tmp)
        ENDIF

!    3.3: P-matrix making by distance or linear shape function
! notes this one: using distance or linear?
            
         nnzi1 = ia_tmp(nnode1+1)-1
         ALLOCATE (iai1(nnode1+1), jai1(nnzi1), Xintp1(nnzi1))
         
!
        iai1(1:nnode1+1) = ia_tmp(1:nnode1+1)
        jai1(1:nnzi1) = ja_tmp(1:nnzi1)
        Xintp1 = 0.d0
!
        IF(ip_inter.EQ.1) THEN
          CALL P_distance(ndim, jmax, nnode1, nnode2, coord1, coord2, nnzi1, iai1, jai1, Xintp1)
        ELSE 
                
          IF(ndim.EQ.2) THEN
          CALL P_linear_2D(jmax,nnode1,nnode2,coord1,coord2,nnzi1,iai1,jai1,Xintp1) 
          ELSE
          CALL P_linear_3D(jmax,nnode1,nnode2,coord1,coord2,nnzi1,iai1,jai1,Xintp1) 
          ENDIF
             
        ENDIF


! 3.4: remove small value
!        IF(MINVAL(Xintp1).LT.alpha) THEN
            
            CALL reduce_CSR_matrix(jmax, nnode1, nnzi1, iai1, jai1, Xintp1, alpha, nnzt)
!
            IF(nnzi1.EQ.nnzt) GOTO 100
!
            nnzi1 = nnzt
            ALLOCATE ( aut(nnzt))

            ja_tmp(1:nnzt) = jai1(1:nnzt)
            aut(1:nnzt) = Xintp1(1:nnzt)
            DEALLOCATE ( jai1, Xintp1)
            ALLOCATE ( jai1(nnzt), Xintp1(nnzt))
            jai1(1:nnzt) = ja_tmp(1:nnzt)
            Xintp1 = aut
            DEALLOCATE ( aut)
!            DEALLOCATE(ia_tmp, ja_tmp)
! 
100 CONTINUE
!
            DEALLOCATE(ia_tmp, ja_tmp)
!        ENDIF
!
!  4:  restriction operator: R = (P)T

            ALLOCATE (iar1(nnode2+1), jar1(nnzi1))
            ALLOCATE (Xrest1(nnzi1))

            CALL mt_trans(nnode1, nnode2, nnzi1, iai1, jai1, iar1, jar1, Xintp1, Xrest1)

!  5:  FOR AC

!      5.1: Connectivity
            CALL connectivity_coarse_nnz(nnode1, nnode2, nnz1, nnzi1, ia1, ja1, &
                                     iai1, jai1, iar1, jar1, nnz2,jmax1)

            ALLOCATE(ia2(nnode2+1), ja2(nnz2), ju2(nnode2))
            
            ia2 = 0
            ja2 = 0
            ju2 = 0
            
!      5.2:  CSR format for coarse grid
!
            CALL connectivity_coarse_CSR(jmax1, nnode1, nnode2, nnz1, nnzi1, ia1, ja1, &
                                     iai1, jai1, iar1, jar1, nnz2, ia2, ja2, ju2)

! - - - -storing 
!ALLOCATE: 
         IF(ilv.EQ.1) THEN
			
		  ALLOCATE(iaiG(nnode1+1), jaiG(nnzi1), iarG(nnode2+1), jarG(nnzi1))
		  ALLOCATE(XintpG(nnzi1), XrestG(nnzi1))
! 
          ALLOCATE(nnodeGC(nlv_glomax),nnzGC(nlv_glomax))
		  ALLOCATE(iaGC(nnode2+1, nlv_glomax), jaGC(nnz2, nlv_glomax), juGC(nnode2, nlv_glomax), auGC(nnz2, nlv_glomax))
		  
		  nnodeGC = 0
		  nnzGC = 0
		  iaGC = 0
		  jaGC = 0
		  juGC = 0
		  auGC = 0.d0
          
!
          ALLOCATE(nnziGC(nlv_glomax))
          nnziGC = 1  ! initial
!
         ALLOCATE(eGt(nnode1),rGt(nnode1))
		  eGt = 0.d0
		  rGt = 0.d0
		  
		 ENDIF
		  
		  IF(ilv.EQ.2) THEN
		  ALLOCATE(iaiGC(nnode1+1, nlv_glomax), jaiGC(nnzi1, nlv_glomax), iarGC(nnode2+1, nlv_glomax), jarGC(nnzi1, nlv_glomax))
		  ALLOCATE(XintpGC(nnzi1, nlv_glomax), XrestGC(nnzi1, nlv_glomax))
		  
		  iaiGC = 0
		  jaiGC = 0
		  iarGC = 0
		  jarGC = 0
		  XintpGC = 0.d0
		  XrestGC = 0.d0
!		  
		  ALLOCATE(eGC(nnode1+1, nlv_glomax),rGC(nnode1+1, nlv_glomax))
		  eGC = 0.d0
		  rGC = 0.d0
		  
          ENDIF		  
		  
! 6: adding to array
			
			IF(ilv.EQ.1) THEN
			nnziG = nnzi1
			iaiG = iai1
			jaiG = jai1
			iarG = iar1
			jarG = jar1	
            XintpG = Xintp1
			XrestG = Xrest1	
			
			ELSE 
			
			nnziGC(ilv) = nnzi1
			iaiGC(1:nnode1+1, ilv) = iai1(1:nnode1+1)
			jaiGC(1:nnzi1, ilv) = jai1(1:nnzi1)
			iarGC(1:nnode2+1, ilv) = iar1(1:nnode2+1)
			jarGC(1:nnzi1, ilv) = jar1(1:nnzi1)	
            XintpGC(1:nnzi1, ilv) = Xintp1(1:nnzi1)
			XrestGC(1:nnzi1, ilv) =	Xrest1(1:nnzi1)			

			ENDIF
			
			nnodeGC(ilv) = nnode2
			nnzGC(ilv) = nnz2
			iaGC(1:nnode2+1, ilv) = ia2(1:nnode2+1)
			jaGC(1:nnz2, ilv) = ja2(1:nnz2)	
			juGC(1:nnode2, ilv) = ju2(1:nnode2)			
! - - - - - 
            DEALLOCATE (iai1, jai1, Xintp1)
            DEALLOCATE (iar1, jar1, Xrest1)
            DEALLOCATE (imap, icoarse)
!
            nnode1 = nnode2
            nnz1 = nnz2
            DEALLOCATE (coord1)
            ALLOCATE (coord1(ndim, nnode1))
            coord1 = coord2
            DEALLOCATE (coord2)

!   CSR matrix
            DEALLOCATE (ia1, ja1)
            ALLOCATE (ia1(nnode2+1), ja1(nnz2))
            ia1 = ia2
            ja1 = ja2

            DEALLOCATE (ia2, ja2, ju2)

! optimize nlevel:
!            IF ((ioplv .EQ. 1) .AND. (ilv .GE. 3)) THEN
!               CALL opt_level(ilv, nlevel, ialv, i)

!               IF (i .EQ. 1) GOTO 20
!
!            END IF
            
            IF(nnode2.EQ.1) THEN
                nlv_glo = ilv
                GOTO 20
            ENDIF
!
         END DO

20       CONTINUE

        IF(nlv_glo.EQ.1) THEN
		  ALLOCATE(iaiGC(1, 1), jaiGC(1, 1), iarGC(1, 1), jarGC(1, 1))
		  ALLOCATE(XintpGC(1, 1), XrestGC(1, 1))	
          iaiGC = 0
		  jaiGC = 0
		  iarGC = 0
		  jarGC = 0
		  XintpGC = 0.d0
		  XrestGC = 0.d0
          
          ALLOCATE(eGC(1, 1),rGC(1, 1))
		  eGC = 0.d0
		  rGC = 0.d0
        ENDIF	
        
! out put for alloate

      nnodecm = MAXVAL( nnodeGC)
          nnzcm = MAXVAL(nnzGC)
          nnzicm = MAXVAL(nnziGC)

! Deallocate
         DEALLOCATE (ia1, ja1)
         DEALLOCATE (coord1)
! -----------------*---------------------*---!
         RETURN
      END
