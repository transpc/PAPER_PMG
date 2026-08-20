      SUBROUTINE reorder_ilup
!
      USE Zmpi         , ONLY: maxmt,maxmt_lu0,maxmt_lu1,          &
                               ia_a,ja_a,ju_a,iend,                &
                               diag_lu,alu0,alu1,ia0,ia1,ja0,ja1,  &
                               lev_typet,levt,lev_typedt,nnp,      &
                               maxmt2,maxmt_r,au_r,ia_r,ja_r,ju_r, &
                               perm_r,permi_r,index_r
      USE Zzone        , ONLY: ncell_fluid
      USE Zbicg        , ONLY: lev_type,levmpi_type,lev,levmpi
      USE Zcore        , ONLY: myrank,np
      USE Zgradoption  , ONLY: non_orth
      IMPLICIT NONE
!.....Local variables
      INTEGER :: i,ip
      INTEGER :: jj,l,len,isize
!.....Local arrays
      INTEGER :: iend_r(ncell_fluid)
      INTEGER :: iar(ncell_fluid+1)
      INTEGER,DIMENSION(:),ALLOCATABLE :: jar
      INTEGER,DIMENSION(:),ALLOCATABLE :: ia_w,ja_w,idiag_w,lev_w
!
      If(np.eq.1) THEN
         levt=lev
         lev_typet=lev_type
         nnp=lev_typet/10
         lev_typedt=mod(lev_typet,10)
         IF(lev_typedt.lt.0 .or. lev_typedt.gt.2) THEN
            WRITE(*,*) ' Wrong lev_type in somaFlow.in ',lev_typet
            STOP
         ENDIF
!        IF(levt.eq.0) lev_typet=0
      ELSE
         levt=levmpi
         lev_typet=levmpi_type
         IF(lev_typet.lt.0 .or. lev_typet.gt.2) THEN
            WRITE(*,*) ' Wrong lev_typet,levmpi_type in somaFlow.in ',lev_typet,levmpi_type
            STOP
         ENDIF
      ENDIF
      IF(levt.lt.-1) THEN
         WRITE(*,*) ' Wrong lev,levmpi in somaFlow.in ',lev,levmpi
         STOP
      ENDIF
!
      If(np.eq.1) THEN
!........Check for tridiagonal
            i=1
            IF(ia_a(i+1)-ia_a(i).ne.2) GOTO 101
            IF(ja_a(1).ne.i  ) GOTO 101
            IF(ja_a(2).ne.i+1) GOTO 101
         DO i=2,ncell_fluid-1
            IF(ia_a(i+1)-ia_a(i).ne.3) GOTO 101
            jj=ia_a(i)
            IF(ja_a(jj  ).ne.i-1) GOTO 101
            IF(ja_a(jj+1).ne.i  ) GOTO 101
            IF(ja_a(jj+2).ne.i+1) GOTO 101
         ENDDO
            i=ncell_fluid
            IF(ia_a(i+1)-ia_a(i).ne.2) GOTO 101
            jj=ia_a(i)
            IF(ja_a(jj  ).ne.i-1) GOTO 101
            IF(ja_a(jj+1).ne.i  ) GOTO 101
         levt=-2
         nnp=0
         lev_typedt=0
!        write(*,*) 'tridiagonal'
         GOTO 102
101      CONTINUE
!        write(*,*) 'not tridiagonal',i
102      CONTINUE
      ENDIF
!      
      IF(myrank.eq.0) THEN
         IF(levt.eq.-2)         WRITE(*,*) '          Requested for fluid: tridiagonal direct solver.'
         IF(levt.eq.-1)         WRITE(*,*) '          Requested for fluid: direct solver.'
         IF(levt.eq. 0)         WRITE(*,*) '          Requested for fluid: ilu0 preconditionner.'
         IF(levt.eq. 1)         WRITE(*,*) '          Requested for fluid: ilu1 preconditionner.'
         IF(levt.eq. 2)         WRITE(*,*) '          Requested for fluid: ilu2 preconditionner.'
         IF(nnp.ne.0) THEN
                                WRITE(*,*) '          Requested for fluid: Domain level decomposition: ',nnp
            IF(lev_typedt.eq.0) WRITE(*,*) '          Requested for fluid: Domain no reordering.'
            IF(lev_typedt.eq.1) WRITE(*,*) '          Requested for fluid: Domain METIS reordering.'
            IF(lev_typedt.eq.2) WRITE(*,*) '          Requested for fluid: Domain CUTHILL reordering.'
         ELSE
           IF(levt.ne.-2) THEN
            IF(lev_typet.eq.0)  WRITE(*,*) '          Requested for fluid: no reordering.'
            IF(lev_typet.eq.1)  WRITE(*,*) '          Requested for fluid: METIS reordering.'
            IF(lev_typet.eq.2)  WRITE(*,*) '          Requested for fluid: CUTHILL reordering.'
           ENDIF
         ENDIF
      ENDIF
!
!
!.....Must always allocate used to call csr_cg_solver
!
      ALLOCATE(perm_r(ncell_fluid),permi_r(ncell_fluid))
      IF(levt.eq.0) THEN
         IF(lev_typet.eq.0) THEN
!
!...........No reordering, lu0 default
!
            maxmt_lu0=0
            maxmt_lu1=0
            DO i=1,ncell_fluid
               l=ju_a(i)-ia_a(i)
               maxmt_lu0=maxmt_lu0+l
               l=iend(i)-ju_a(i)
               maxmt_lu1=maxmt_lu1+l
            ENDDO
         ELSE
!
!...........Build a till iend without diagonal
!
            iar(1)=1
            DO i=1,ncell_fluid
               l=iend(i)-ia_a(i)
               iar(i+1)=iar(i)+l
            ENDDO
            len=iar(ncell_fluid+1)-1
            ALLOCATE(jar(len))
            maxmt2=len+ncell_fluid
            ip=1
            DO i=1,ncell_fluid
               DO jj=ia_a(i),ju_a(i)-1
                  jar(ip)=ja_a(jj)
                  ip=ip+1
               ENDDO
               DO jj=ju_a(i)+1,iend(i)
                  jar(ip)=ja_a(jj)
                  ip=ip+1
               ENDDO
            ENDDO
!
!...........Choose reordering scheme
!
!           write(*,*) 'lev_typet',lev_typet,lev_typedt,nnp
            IF(nnp.eq.0) THEN
               IF(lev_typet.eq.1) THEN
                  CALL reorder_metis(ncell_fluid,iar,jar, &
                                     perm_r,permi_r,len)
               ELSEIF(lev_typet.eq.2) THEN
                  CALL reorder_cuthill(ncell_fluid,iar,jar, &
                                       perm_r,permi_r,len)
               ENDIF
!              write(*,*) '9==>',maxmt,maxmt2
            ENDIF
            ALLOCATE(ia_r(ncell_fluid+1))
            ALLOCATE(ju_r(ncell_fluid))
            ALLOCATE(ja_r(maxmt2))
            ALLOCATE(index_r(maxmt2))
            CALL final_r(ncell_fluid,maxmt,maxmt2, &
                         ia_a,iend,ja_a,            &
                         ia_r,ju_r,ja_r,            &
                         perm_r,permi_r,index_r)
            maxmt_lu0=0
            maxmt_lu1=0
            DO i=1,ncell_fluid
               l=ju_r(i)-ia_r(i)
               maxmt_lu0=maxmt_lu0+l
               l=ia_r(i+1)-1-ju_r(i)
               maxmt_lu1=maxmt_lu1+l
            ENDDO
!              write(*,*) '9==>',maxmt_lu0,maxmt_lu1
         ENDIF
      ELSEIF(levt.ge.-1) THEN
         IF(lev_typet.eq.0) THEN
!
!...........No reordering, lup ia_r,ju_r,ja_r array till iend
!
            ALLOCATE(ia_r(ncell_fluid+1))
            ia_r(1)=1
            DO i=1,ncell_fluid
               l=iend(i)-ia_a(i)+1
               ia_r(i+1)=ia_r(i)+l
            ENDDO
            maxmt2=ia_r(ncell_fluid+1)-1
            ALLOCATE(ju_r(ncell_fluid))
            ALLOCATE(ja_r(maxmt2))
            ip=1
            DO i=1,ncell_fluid
               DO jj=ia_a(i),iend(i)
                  ja_r(ip)=ja_a(jj)
                  IF(ja_a(jj).eq.i) THEN
                     ju_r(i)=ip
                     ip=ip+1
                     EXIT
                  ENDIF
                  ip=ip+1
               ENDDO
               DO jj=ju_a(i)+1,iend(i)
                  ja_r(ip)=ja_a(jj)
                  ip=ip+1
               ENDDO
            ENDDO
         ELSE
!
!...........Build a till iend without diagonal
!
            iar(1)=1
            DO i=1,ncell_fluid
               l=iend(i)-ia_a(i)
               iar(i+1)=iar(i)+l
            ENDDO
            len=iar(ncell_fluid+1)-1
            ALLOCATE(jar(len))
            maxmt2=len+ncell_fluid
            ip=1
            DO i=1,ncell_fluid
               DO jj=ia_a(i),ju_a(i)-1
                  jar(ip)=ja_a(jj)
                  ip=ip+1
               ENDDO
               DO jj=ju_a(i)+1,iend(i)
                  jar(ip)=ja_a(jj)
                  ip=ip+1
               ENDDO
            ENDDO
!
!...........Choose reordering scheme
!
!           write(*,*) 'lev_typet',lev_typet,lev_typedt,nnp
            IF(nnp.eq.0) THEN
               IF(lev_typet.eq.1) THEN
                  CALL reorder_metis(ncell_fluid,iar,jar, &
                                     perm_r,permi_r,len)
               ELSEIF(lev_typet.eq.2) THEN
                  CALL reorder_cuthill(ncell_fluid,iar,jar, &
                                       perm_r,permi_r,len)
               ENDIF
!              write(*,*) '9==>',maxmt,maxmt2
            ELSE
            ENDIF
            ALLOCATE(ia_r(ncell_fluid+1))
            ALLOCATE(ju_r(ncell_fluid))
            ALLOCATE(ja_r(maxmt2))
            ALLOCATE(index_r(maxmt2))
            CALL final_r(ncell_fluid,maxmt,maxmt2, &
                         ia_a,iend,ja_a,            &
                         ia_r,ju_r,ja_r,            &
                         perm_r,permi_r,index_r)
!
         ENDIF
!
!........Do symbolic factorisation in ia_w,idiag_w,ja_w
!
!        write(*,*) '==>',isize,maxmt2
         isize=200*maxmt2
!        isize = 2000000000 
         write(*,*) '0-==>',isize,maxmt2
         ALLOCATE(lev_w(isize))
         ALLOCATE(ia_w(ncell_fluid+1))
         ALLOCATE(ja_w(isize))
         ALLOCATE(idiag_w(ncell_fluid))
         IF(levt.eq.-1) THEN
            CALL dgefs2s(ncell_fluid,maxmt2,isize, &
                         ia_r,ja_r,                &
                         ia_w,ja_w,idiag_w)
         ELSE
            CALL dgefs2s_lev(ncell_fluid,maxmt2,isize,     &
                             ia_r,ja_r,                    &
                             ia_w,ja_w,idiag_w,lev_w,levt)
         ENDIF
         IF(myrank.eq.0) THEN
!           WRITE(*,*) 'fill in ratio : ',ia_w(ncell_fluid+1)-1,ia_r(ncell_fluid+1)-1
            WRITE(*,401) 'fill in ratio : ',dble(ia_w(ncell_fluid+1)-1)/(ia_r(ncell_fluid+1)-1)
         ENDIF
!
!........Size with fill-in computed put back in ia_r,ju_r,ha_r
!
         maxmt_r=ia_w(ncell_fluid+1)-1
         DEALLOCATE(ja_r)
         ALLOCATE(ja_r(maxmt_r))
         ALLOCATE(au_r(maxmt_r))
         do i=1,ncell_fluid
            ju_r(i)=idiag_w(i)
            iend_r(i)=ia_w(i+1)-1
         enddo
         do i=1,ncell_fluid+1
            ia_r(i)=ia_w(i)
         enddo
         do i=1,maxmt_r
            ja_r(i)=ja_w(i)
         enddo
         maxmt_lu0=0
         maxmt_lu1=0
         DO i=1,ncell_fluid
            l=ju_r(i)-ia_r(i)
            maxmt_lu0=maxmt_lu0+l
            l=iend_r(i)-ju_r(i)
            maxmt_lu1=maxmt_lu1+l
         ENDDO
      ELSEIF(levt.eq.-2) THEN
          IF(non_orth.gt.0)THEN
             ALLOCATE(diag_lu(ncell_fluid))
             ALLOCATE(alu0(ncell_fluid),alu1(ncell_fluid))
          ENDIF
          RETURN
      ENDIF
401   FORMAT(11x,a,f6.2,2(i5,2x))
!
      ALLOCATE(diag_lu(ncell_fluid))
      ALLOCATE(alu0(maxmt_lu0),alu1(maxmt_lu1))
      ALLOCATE(ja0(maxmt_lu0),ja1(maxmt_lu1))
      ALLOCATE(ia0(ncell_fluid+1),ia1(ncell_fluid+1))
!
      END SUBROUTINE reorder_ilup
