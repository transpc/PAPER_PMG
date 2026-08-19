      SUBROUTINE PREP_GMG

         USE MD_parameter, ONLY: nv_max, ndim, teta,teta_p, nf_max, alpha, ip_inter
         USE MD_geometry, ONLY: nelem, coord, num_neigh_mg, neigh_mg, imap
         USE MD_MG_coord
         USE MD_MG_index, ONLY: nlevel, mxnbne,ip_nmax, nmax, AR_hi, ioplv,ip_lev, report_text,nlevel_N
         USE MD_MG_matrix
         USE MD_matrix, ONLY: nnz, ia, ja, ju, au
         USE MD_connectivity, ONLY: ia_neigh, ja_neigh, nnz_neigh,   &
                                    ia_neighc, ja_neighc, nnz_neighc, nnz_tmp, ia_tmp, ja_tmp
         USE MD_MG_Global_C, ONLY: nlv_glo
! ---
         IMPLICIT NONE

! ---temp
         INTEGER(4) i, k, nv, nd, j, k1, k2, l, i1, i2, nnd,j1
         INTEGER(4) ilv, nnzt, ntmp, ncolc1, ncolc2
         REAL(8) xtmp1(3), xtmp2(3), dist, dmin
         INTEGER(4) jmax, jmax1
         INTEGER(4):: alstatus
         REAL(8) teta1!alpha,

!         INTEGER(4), DIMENSION(:), ALLOCATABLE :: nneit!,imapt
!         INTEGER(4), DIMENSION(:, :), ALLOCATABLE :: ineit
!         INTEGER(4), DIMENSION(:), ALLOCATABLE :: iwk      !, iat, jat
!         INTEGER(4), DIMENSION(:, :), ALLOCATABLE :: iwork
!         REAL(8), DIMENSION(:), ALLOCATABLE :: aut

         CHARACTER(LEN=10) :: filename
!
! initial
         ALLOCATE (ialv(nlevel+1))
         ialv(1) = 1
         ialv(2) = ialv(1)+nelem

!
         nnode1 = nelem
         nnz1 = nnz
         ALLOCATE (coord1(ndim, 1))
         coord1 = 0.d0
         ALLOCATE (ia1(1), ja1(1))
         ia1 = 0
         ja1 = 0

         ncolc1 = 0
         ncolc2 = 0

         nmax1 = nf_max

!         DEALLOCATE (coord)
! initial
         ALLOCATE (iac(1), jac(1))
         ALLOCATE (iai(1), jai(1), Xintp(1))
         ALLOCATE (iar(1), jar(1), Xrest(1))
         iac = 1
         jac = 1
         iai = 1
         jai = 1
         iar = 1
         jar = 1
         Xintp = 1.d0
         Xrest = 1.d0

         ALLOCATE (imapc(1), icoarsef(1), inmax(nlevel))
         imapc = 0
         icoarsef = 0
         inmax(1) = nmax1
! for FVM
! ia_neigh, ja_neigh
         nnz_neigh = SUM(num_neigh_mg(1:nelem))
         ALLOCATE(ia_neigh(nelem+1), ja_neigh(nnz_neigh))
         
         ia_neigh(1) = 1
         
         DO i = 1, nelem
             j = num_neigh_mg(i)
             k = ia_neigh(i)
             DO j1 = 1, j
                 ja_neigh(k) = neigh_mg(j1,i)
                 k = k+1
             ENDDO
             ia_neigh(i+1) = k
         ENDDO
         
         IF(k.NE.(nnz_neigh+1)) THEN
             WRITE(*,*)'error in ia_neigh'
             STOP
             
         ENDIF
         
!

         DO ilv = 2, nlevel

! 1: neighbor nodes of each nodes  - - - - - - - - - - !

!            IF (ilv .GT. 2) THEN
!               ALLOCATE (nnei(nnode1), inei(nmax1, nnode1))
!
!               nnei(1:nnode1) = nneit(1:nnode1)
!               inei(1:nmax1, 1:nnode1) = ineit(1:nmax1, 1:nnode1)
!               DEALLOCATE (nneit, ineit)

!            END IF

! 2: coarsening step
            teta1 = teta
            IF (AR_hi .EQ. 1) THEN
               IF (ilv .GE. 4) teta1 = teta+0.1
               IF (ilv .GE. 6) teta1 = teta+0.2
               IF (ilv .GE. 9) teta1 = teta+0.3
               IF (teta1 .GT. 0.92) teta1 = 0.92
            END IF

            nnode2 = INT(0.6*nnode1)
            IF (nnode1 .LE. 2000) nnode2 = nnode1

            ALLOCATE (imap(nnode2), icoarse(nnode1),stat=alstatus)
            
            IF (alstatus/=0) THEN
               report_text = 'not enough memory,PREP_GMG_FEM, imap'
               CALL STOP_MPI(report_text)
            ENDIF
            
            imap = 0
            icoarse = 0
            
            IF (ilv .EQ. 2) THEN

               CALL coarsening_semi(ndim, nnode1, nmax1, nnz_neigh, ia_neigh, ja_neigh, nnode2, imap, icoarse, teta, coord)

            ELSE
               CALL coarsening_semi(ndim, nnode1, nmax1, nnz_neighc, ia_neighc, ja_neighc, nnode2, imap, icoarse, teta1, coord1)

            END IF
!
            ALLOCATE (coord2(ndim, nnode2))
!
            IF(ilv.EQ.2) THEN
                
              DO i = 1, nnode2
                 j = imap(i)

                 coord2(1:ndim, i) = coord(1:ndim, j)

              END DO
            ELSE
            DO i = 1, nnode2
               j = imap(i)

               coord2(1:ndim, i) = coord1(1:ndim, j)

            END DO
            ENDIF
            

! 3: Iterpolation procedure
!            jmax = nmax1
!    3.1: iwk-neighbor nodes of each fine-cell
! finding nnz_neighc            
            IF (ilv .EQ. 2) THEN
              IF(ip_lev.EQ.1) THEN
               CALL neighbor_fine_graph_nnz(nnode1,nnode2,nnz_neigh, ia_neigh, ja_neigh,icoarse,nnz_tmp, jmax)
              ELSE
               CALL neighbor_fine_graph2_nnz(nnode1,nnode2,nnz_neigh, ia_neigh, ja_neigh,icoarse,nnz_tmp, jmax)
              ENDIF
              
            ELSE
              IF(ip_lev.EQ.1) THEN 
               CALL neighbor_fine_graph_nnz(nnode1,nnode2,nnz_neighc, ia_neighc, ja_neighc,icoarse,nnz_tmp, jmax)
              ELSE
               CALL neighbor_fine_graph2_nnz(nnode1,nnode2,nnz_neighc, ia_neighc, ja_neighc,icoarse,nnz_tmp, jmax)
              ENDIF
              
            END IF
            
! allocate:
!                         write(999,*) 'test1', alstatus
!            IF(ALLOCATED(ia_tmp))      DEALLOCATE(ia_tmp)
!            IF(ALLOCATED(ja_tmp))      DEALLOCATE(ja_tmp)


            ALLOCATE(ia_tmp(nnode1+1), ja_tmp(nnz_tmp),stat=alstatus)

            IF (alstatus/=0) THEN
               report_text = 'not enough memory,PREP_GMG_FEM, ja_tmp'
             write(999,*) nnz_tmp,nnode1,jmax
               CALL STOP_MPI(report_text)
            ENDIF

            ia_tmp = 0
            ja_tmp = 0
! fiding neighors: ja_tmp
            
!           
            IF (ilv .EQ. 2) THEN
              IF(ip_lev.EQ.1) THEN
               CALL neighbor_fine_graph(nnode1,nnode2,nnz_neigh, ia_neigh, ja_neigh,icoarse,nnz_tmp, jmax, ia_tmp, ja_tmp)
              ELSE
               CALL neighbor_fine_graph2(nnode1,nnode2,nnz_neigh, ia_neigh, ja_neigh,icoarse,nnz_tmp, jmax, ia_tmp, ja_tmp)
              ENDIF
              
            DEALLOCATE(ia_neigh, ja_neigh)
            
            ELSE
              IF(ip_lev.EQ.1) THEN 
               CALL neighbor_fine_graph(nnode1,nnode2,nnz_neighc, ia_neighc, ja_neighc,icoarse,nnz_tmp, jmax, ia_tmp, ja_tmp)
              ELSE
               CALL neighbor_fine_graph2(nnode1,nnode2,nnz_neighc, ia_neighc, ja_neighc,icoarse,nnz_tmp, jmax, ia_tmp, ja_tmp)
              ENDIF
              
            DEALLOCATE(ia_neighc, ja_neighc)
            END IF

!    3.2: reduce to n-max nodes by order of distance

!    3.2: reduce iwk (using teta)

        i = 0
        IF(teta_p.GT.0.1) i = i + 1 
        IF((ip_nmax.NE.0).AND.(jmax.GT.ip_nmax)) i = i + 1     
!
        IF(i.NE.0) THEN
        IF (ilv .EQ. 2) THEN
        CALL reduce_neibor(ndim,jmax,ip_nmax,teta_p,nnode1,nnode2,coord,coord2,nnz_tmp,ia_tmp,ja_tmp)            
        ELSE    
        CALL reduce_neibor(ndim,jmax,ip_nmax,teta_p,nnode1,nnode2,coord1,coord2,nnz_tmp,ia_tmp,ja_tmp)
        ENDIF
        ENDIF
!    3.3: P-matrix making by distance or linear shape function

        nnzi1 = ia_tmp(nnode1+1)-1
        ALLOCATE (iai1(nnode1+1), jai1(nnzi1), Xintp1(nnzi1),stat=alstatus )

        IF (alstatus/=0) THEN
         report_text = 'not enough memory,PREP_GMG_FEM, Xintp1'
         CALL STOP_MPI(report_text)
        ENDIF
        
        iai1(1:nnode1+1) = ia_tmp(1:nnode1+1)
        jai1(1:nnzi1) = ja_tmp(1:nnzi1)
        Xintp1 = 0.d0 
        
! - - - - - - - - - - - - - - - 
      IF(ilv.EQ.2) THEN
        IF(ip_inter.EQ.1) THEN
          CALL P_distance(ndim, jmax, nnode1, nnode2, coord, coord2, nnzi1, iai1, jai1, Xintp1)
        ELSE 
                
          IF(ndim.EQ.2) THEN
          CALL P_linear_2D(jmax,nnode1,nnode2,coord,coord2,nnzi1,iai1,jai1,Xintp1) 
          ELSE
          CALL P_linear_3D(jmax,nnode1,nnode2,coord,coord2,nnzi1,iai1,jai1,Xintp1) 
          ENDIF
             
        ENDIF
        
      ELSE
            
        IF(ip_inter.EQ.1) THEN
          CALL P_distance(ndim, jmax, nnode1, nnode2, coord1, coord2, nnzi1, iai1, jai1, Xintp1)
        ELSE 
                
          IF(ndim.EQ.2) THEN
          CALL P_linear_2D(jmax,nnode1,nnode2,coord1,coord2,nnzi1,iai1,jai1,Xintp1) 
          ELSE
          CALL P_linear_3D(jmax,nnode1,nnode2,coord1,coord2,nnzi1,iai1,jai1,Xintp1) 
          ENDIF
             
        ENDIF
        
      ENDIF

! 3.4: remove small value
!      alpha = 0.005
            CALL reduce_CSR_matrix(jmax,nnode1, nnzi1, iai1, jai1, Xintp1, alpha, nnzt)
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
! 
100 CONTINUE

            DEALLOCATE(ia_tmp, ja_tmp)

!  4:  restriction operator: R = (P)T

            ALLOCATE (iar1(nnode2+1), jar1(nnzi1))
            ALLOCATE (Xrest1(nnzi1))

            CALL mt_trans(nnode1, nnode2, nnzi1, iai1, jai1, iar1, jar1, Xintp1, Xrest1)

!  5:  FOR AC

!      5.1: Connectivity-nnz

            IF(ilv.EQ.2) THEN
            CALL connectivity_coarse_nnz(nnode1, nnode2, nnz1, nnzi1, ia, ja, &
                                     iai1, jai1, iar1, jar1, nnz2,jmax1)                
            ELSE
            CALL connectivity_coarse_nnz(nnode1, nnode2, nnz1, nnzi1, ia1, ja1, &
                                     iai1, jai1, iar1, jar1, nnz2,jmax1)                
            ENDIF
!
            
            ALLOCATE(ia2(nnode2+1), ja2(nnz2), ju2(nnode2),stat=alstatus)
            
            IF (alstatus/=0) THEN
            report_text = 'not enough memory,PREP_GMG_FEM, ia2,ja2'
            CALL STOP_MPI(report_text)
            ENDIF
            ia2 = 0
            ja2 = 0
            ju2 = 0
!
!      5.2:  CSR format for coarse grid
          IF(ilv.EQ.2) THEN
            CALL connectivity_coarse_CSR(jmax1, nnode1, nnode2, nnz1, nnzi1, ia, ja, &
                                     iai1, jai1, iar1, jar1, nnz2, ia2, ja2, ju2)
          ELSE
            CALL connectivity_coarse_CSR(jmax1, nnode1, nnode2, nnz1, nnzi1, ia1, ja1, &
                                     iai1, jai1, iar1, jar1, nnz2, ia2, ja2, ju2)
          ENDIF


! 6: adding to array
            ialv(ilv+1) = ialv(ilv)+nnode2

! for Ac
            ntmp = ncolc2+nnode2
            nnzt = iac(ncolc2+1)-1
            ALLOCATE (iat(ntmp+1), jat(nnzt+nnz2))
            iat(1:ncolc2+1) = iac(1:ncolc2+1)
            iat(ncolc2+2:ntmp+1) = ia2(2:nnode2+1)+nnzt
            jat(1:nnzt) = jac(1:nnzt)
            jat(nnzt+1:nnzt+nnz2) = ja2(1:nnz2)+ncolc2
!      jut (1:ncolc2) = juc (1:ncolc2)
!      jut (ncolc2+1:ncolc2+nnode2) = ju2 (1:nnode2) + nnzt

            DEALLOCATE (iac, jac)

            nnzt = nnzt+nnz2
            ALLOCATE (iac(ntmp+1), jac(nnzt))
            iac = iat
            jac = jat
            DEALLOCATE (iat, jat)

! for Xintp
            ntmp = ncolc1+nnode1
            nnzt = iai(ncolc1+1)-1
            ALLOCATE (iat(ntmp+1), jat(nnzt+nnzi1), aut(nnzt+nnzi1))
            iat(1:ncolc1+1) = iai(1:ncolc1+1)
            iat(ncolc1+2:ntmp+1) = iai1(2:nnode1+1)+nnzt
            jat(1:nnzt) = jai(1:nnzt)
            jat(nnzt+1:nnzt+nnzi1) = jai1(1:nnzi1)+ncolc2
            aut(1:nnzt) = Xintp(1:nnzt)
            aut(nnzt+1:nnzt+nnzi1) = Xintp1(1:nnzi1)

            DEALLOCATE (iai, jai, Xintp)

            nnzt = nnzt+nnzi1
            ALLOCATE (iai(ntmp+1), jai(nnzt), Xintp(nnzt))
            iai = iat
            jai = jat
            Xintp = aut
            DEALLOCATE (iat, jat, aut)
            DEALLOCATE (iai1, jai1, Xintp1)

! for Xrest
            ntmp = ncolc2+nnode2
            nnzt = iar(ncolc2+1)-1
            ALLOCATE (iat(ntmp+1), jat(nnzt+nnzi1), aut(nnzt+nnzi1))
            iat(1:ncolc2+1) = iar(1:ncolc2+1)
            iat(ncolc2+2:ntmp+1) = iar1(2:nnode2+1)+nnzt
            jat(1:nnzt) = jar(1:nnzt)
            jat(nnzt+1:nnzt+nnzi1) = jar1(1:nnzi1)+ncolc1
            aut(1:nnzt) = Xrest(1:nnzt)
            aut(nnzt+1:nnzt+nnzi1) = Xrest1(1:nnzi1)

            DEALLOCATE (iar, jar, Xrest)

            nnzt = nnzt+nnzi1
            ALLOCATE (iar(ntmp+1), jar(nnzt), Xrest(nnzt))
            iar = iat
            jar = jat
            Xrest = aut
            DEALLOCATE (iat, jat, aut)
            DEALLOCATE (iar1, jar1, Xrest1)

! imap,
            ntmp = ncolc2+nnode2
            ALLOCATE (imapt(ntmp))
            imapt(1:ncolc2) = imapc(1:ncolc2)
            imapt(ncolc2+1:ntmp) = imap(1:nnode2)

            DEALLOCATE (imapc)

            ALLOCATE (imapc(ntmp))
            imapc = imapt

            DEALLOCATE (imapt, imap)

! icoarse
            ntmp = ncolc1+nnode1
            ALLOCATE (imapt(ntmp))
            imapt(1:ncolc1) = icoarsef(1:ncolc1)
            imapt(ncolc1+1:ntmp) = icoarse(1:nnode1)

            DEALLOCATE (icoarsef)

            ALLOCATE (icoarsef(ntmp))
            icoarsef = imapt

            DEALLOCATE (imapt, icoarse)

! the coarest level:
!      IF(ilv .EQ. nlevel) THEN
!          ALLOCATE(aus(nnz2),alus(nnz2))
!          ALLOCATE(ias(nnode2+1),jas(nnz2),jus(nnode2))
!          ias = ia2
!          jas = ja2
!          jus = ju2
!      ENDIF

! 7: Updating for next level:
            ncolc1 = ncolc1+nnode1
            ncolc2 = ncolc2+nnode2

            nnode1 = nnode2
            nnz1 = nnz2
            DEALLOCATE (coord1)
            ALLOCATE (coord1(ndim, nnode1))
            coord1 = coord2
            DEALLOCATE (coord2)

!   neighbor nodes

            nmax2 = 0

            DO i = 1, nnode2
               i1 = ia2(i)
               i2 = ia2(i+1)-1

               nnd = i2-i1
               nmax2 = MAX(nmax2, nnd)
            END DO

            IF (nmax2 .EQ. 0) nmax2 = 1
            
            inmax(ilv) = nmax2
            nmax1 = nmax2

! neighbor for next level:
            nnz_neighc = nnz2 - nnode2                         ! only for neighbor nodes
            
            ALLOCATE (ia_neighc(nnode2+1), ja_neighc(nnz_neighc))

            ia_neighc(1) = 1
            nnzt = 0
            DO i = 1, nnode2
               i1 = ia2(i)
               i2 = ia2(i+1)-1

               nnd = 0
               DO j = i1, i2
                  k = ja2(j)
                  IF (k .EQ. i) CYCLE
                  nnd = nnd+1
                  nnzt = nnzt+1
                  ja_neighc(nnzt) = k
               END DO
               ia_neighc(i+1) = ia_neighc(i)+ nnd

            END DO



!   CSR matrix
            DEALLOCATE (ia1, ja1)
            ALLOCATE (ia1(nnode2+1), ja1(nnz2))
            ia1 = ia2
            ja1 = ja2

            DEALLOCATE (ia2, ja2, ju2)

! optimize nlevel:
            
! only 1 node remaining
            
            IF(nnode2.EQ.1) THEN
               CALL opt_level(ilv, nlevel, ialv, i)
               GOTO 20
            ENDIF
            
            IF ((ioplv .EQ. 1) .AND. (ilv .GE. 3)) THEN
               CALL opt_level(ilv, nlevel, ialv, i)
               
               nlevel_N = nlevel_N-nlv_glo

               IF (i .EQ. 1) GOTO 20

            END IF

!
         END DO

20       CONTINUE

!
!     nnods = ialv(nlevel+1)-ialv(nlevel)
!     ncolf = ialv(nlevel+1)-ialv(1)
!     ncolc = ncolf - nelem

!     ALLOCATE(r(nelem),rt(ncolf),rc(ncolc))
!     ALLOCATE(rs(nnods),es(nnods),e(ncolc),et(ncolf))

!     nnzt = iac(ncolc+1)-1
!     ALLOCATE(auc(nnzt))
!     auc = 0.d0

! Deallocate
         DEALLOCATE (coord)
         DEALLOCATE (ia_neighc, ja_neighc)
         DEALLOCATE (ia1, ja1)
         DEALLOCATE (coord1)
! -----------------*---------------------*---!
         RETURN
      END

! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
! = = = = = = = = = = = = = = = = = = = = = = = = =

! - - - - - - - - - - - - - - - - - - - - - - - !
      !

! = = = = = = = = = = = = = = = = = = = = = = = = = !
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

      SUBROUTINE reduce_CSR_matrix(jmax,ne, nnz, ia, ja, au, alpha, nnz1)

!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
!       remove the small value of restriction operator                 !
!      after that, the matrix is more sparse                           !

         IMPLICIT NONE
!
         INTEGER ne, nnz,jmax
         INTEGER ia(ne+1), ja(nnz)
         REAL*8 alpha
         REAL*8 au(nnz)
! out
         INTEGER nnz1
!
         INTEGER ie, i, i1, i2, nd, ni(jmax)
         REAL*8 xt, ux(jmax)
         INTEGER(4), DIMENSION(:), ALLOCATABLE :: ia1

!
         ALLOCATE (ia1(ne+1))

         ia1 = 0

         ia1(1) = 1
         nnz1 = 0

         DO ie = 1, ne
            i1 = ia(ie)
            i2 = ia(ie+1)-1

            DO i = i1, i2
               xt = au(i)

               IF (abs(xt) .ge. alpha) THEN
                  nnz1 = nnz1+1
                  ja(nnz1) = ja(i)
                  au(nnz1) = xt
               END IF

            END DO

            ia1(ie+1) = nnz1+1

         END DO

! set values
         ia = ia1

         DEALLOCATE (ia1)
         
! New for re-ordering JA, au
         DO ie = 1, ne
            i1 = ia(ie)
            i2 = ia(ie+1)-1

            DO i = i1, i2
               nd = i2-i1+1
               IF(nd.EQ.1) CYCLE
               ni(1:nd) = ja(i1:i2)
               ux(1:nd) = au(i1:i2)
               CALL bubble_sort_2(nd,ni,ux)
! 
               ja(i1:i2) = ni(1:nd)
               au(i1:i2) = ux(1:nd)
            ENDDO
            
         ENDDO   
            
         RETURN
    END

! = = = = = = = = = = = = = = = = = = = = = = = = = = = 
! = = = = = = = = = = = = = = = = = = = = = = = = = = !
      SUBROUTINE opt_level(ilv, nlevel, ialv, id)

         USE MD_parameter, ONLY: ndom
         USE MD_MG_index, ONLY: nlevel_N, n1_min, n2_min

         IMPLICIT NONE

         INTEGER(4) ilv, nlevel
         INTEGER(4) ialv(nlevel+1)
         INTEGER(4) id

!
         INTEGER(4) n0, n1, n2

         id = 0

!         n0 = ialv(ilv)-ialv(ilv-1)
         n1 = ialv(ilv+1)-ialv(ilv)

!         n2 = INT(n1/ndom)
         ! test
!   write(*,*)'n0,',n0,n1,n2,ndom

!         IF ((n2 .LE. n2_min) .OR. (n1 .LE. n1_min)) THEN
         IF (n1 .LE. n1_min) THEN
            id = 1
            nlevel_N = ilv
         END IF

         RETURN

      END

