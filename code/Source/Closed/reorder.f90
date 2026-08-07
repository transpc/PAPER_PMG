!
      SUBROUTINE reorder_cuthill(n,iar,jar,      &
                                 perm,permi,len)
!      
      IMPLICIT NONE
!
!.....Input
      INTEGER :: n,len
      INTEGER :: iar(n+1),jar(len)
!.....Output
      INTEGER :: perm(n),permi(n)
!.....Local variables
!.....Local arrays
!
      CALL cuthill(n,iar,jar,perm,permi)
!     do i=1,n
!        write(*,*) i,perm(i),permi(i)
!     enddo
!     stop 10
!
      END SUBROUTINE reorder_cuthill
!
!-----------------------------------------------------------------------
!
      SUBROUTINE reorder_metis(n,iar,jar,      &
                               perm,permi,len)
!
      IMPLICIT NONE
!.....Input
      INTEGER :: n,len
      INTEGER :: iar(n+1),jar(len)
!.....Output
      INTEGER :: perm(n),permi(n)
!.....Local variables
!.....Local arrays
      INTEGER :: moptions(40)
      INTEGER :: vwgt(n)
!
!DEC$IF defined (metis_flag)
!
!
      vwgt(:)=1
      CALL metis_setdefaultoptions(moptions)
      moptions(17)=1
      CALL METIS_NodeND(n,iar,jar,vwgt,moptions,perm,permi)
!DEC$ELSE
!.....In windows debug no metis use cuthill
      WRITE(*,*) '          WINDOWS debug METIS not available: switch to cuthill'
      CALL cuthill(n,iar,jar,perm,permi)
!     vwgt(:)=1
!     CALL metis_setdefaultoptions(moptions)
!     moptions(17)=1
!     CALL METIS_NodeND(n,iar,jar,vwgt,moptions,perm,permi)
!DEC$ENDIF
!
      END SUBROUTINE reorder_metis
!
!-----------------------------------------------------------------------
!
      SUBROUTINE reorder_domain(nelem,maxmt,i_neigh0,j_neigh0,j_nbcon0,n,maxmt1,ia0,ja0,ju,iend, &
                             perm,permi,len,celem,lev_typedt,nnp,cinterp,ncell_fluid1)
!
!DEC$IF defined (metis_flag)
!      
      IMPLICIT NONE
!
!.....Input
      INTEGER :: n,maxmt1,maxmt
      INTEGER :: ia0(n+1),ju(n),iend(n)
      INTEGER :: ja0(maxmt1)
      INTEGER :: i_neigh0(nelem+1),j_neigh0(maxmt),j_nbcon0(maxmt) 
      INTEGER :: nelem
      INTEGER :: celem(nelem)
      INTEGER :: cinterp(nnp)
      INTEGER :: lev_typedt
!.....Output
      INTEGER :: len
      INTEGER :: perm(nelem),permi(nelem)
!     local variables
      INTEGER :: i,j,k,l,ip,ie,jp,iptr,jj,j0,j1,k1
      INTEGER :: ip1,ip2
      INTEGER :: cinter
      INTEGER :: ncell_fluid
      INTEGER :: cext0
      INTEGER :: nnp
!     PARAMETER(nnp=4)
!     INTEGER :: cinterp(nnp)
!     local arrays
      INTEGER :: celem_offset
      INTEGER :: irecv_cnt(nnp),jsend_cnt(nnp)
      INTEGER :: flag(nelem),flagt(nelem)
      INTEGER :: ia(nnp+1),ja(nelem)
      INTEGER :: icount(nnp)
      INTEGER :: ncell_fluid1(nnp)
      INTEGER :: moptions(40)
      INTEGER :: edgecut(nelem)
      INTEGER :: ncon(nnp),vwgt(nelem*1)
      INTEGER,DIMENSION(:),ALLOCATABLE :: adjwgt
      INTEGER :: jperms(nelem)
      INTEGER :: ia_sub(nelem+1)
      REAL(4) :: tpwgts(nnp),ubvec(1)
!
      INTEGER,ALLOCATABLE :: ja_sub(:)
      INTEGER,DIMENSION(:),ALLOCATABLE :: jperm
!!!!!!!!!!!!!!!
      INTEGER :: maxmt2
      INTEGER,DIMENSION(:),ALLOCATABLE :: ia2,ja2,perm2,permi2
      INTEGER :: niut  
      INTEGER,DIMENSION(:),ALLOCATABLE :: iut,si,jsend
!
!.....Build a till iend without diagonal
      ia_sub(1)=1
      DO i=1,nelem
         l=iend(i)-ia0(i)
         ia_sub(i+1)=ia_sub(i)+l
      ENDDO
      len=ia_sub(nelem+1)-1
      ALLOCATE(ja_sub(len))
      len=len+nelem
      ip=1
      DO i=1,nelem
         DO j=ia0(i),ju(i)-1
            ja_sub(ip)=ja0(j)
            ip=ip+1
         ENDDO
         DO j=ju(i)+1,iend(i)
            ja_sub(ip)=ja0(j)
            ip=ip+1
         ENDDO
      ENDDO
!      
!.....Domain partition by METIS
!
      DO i=1,nelem
         celem(i)=0
      ENDDO
!     IF (metis.eq.1) THEN
         vwgt(:)=1
         ncon(:)=1
         adjwgt(:)=1
         tpwgts(:)=1.0/FLOAT(nnp)
         ubvec(:)=1.001
         CALL metis_setdefaultoptions(moptions)
         moptions(17)=1
         CALL METIS_PartGraphKway(nelem,1,ia_sub,ja_sub,vwgt,0,adjwgt,nnp,tpwgts,ubvec,moptions,edgecut,celem) ! METIS-5.0,
         celem_offset=1-MINVAL(celem(1:nelem))
         celem=celem+celem_offset
!     ENDIF
      DO i=1,nelem
!        write(*,*) i,celem(i)
      ENDDO
      DO i=1,nelem
         perm(i)=0
         permi(i)=0
      ENDDO
!
!     build index based on proc id to access global cell array directly
!     get ncell_fluid for all procs
!
      ia(1)=1
      DO ip=1,nnp
         ia(ip+1)=0
         icount(ip)=0
      ENDDO
      DO ie=1,nelem
         ip=celem(ie)
         ia(ip+1)=ia(ip+1)+1
      ENDDO
      DO ip=1,nnp
         ncell_fluid1(ip)=ia(ip+1)
         ia(ip+1)=ia(ip+1)+ia(ip)
      ENDDO
      DO ie=1,nelem
         ip=celem(ie)
         ja(ia(ip)+icount(ip))=ie
         icount(ip)=icount(ip)+1
      ENDDO
!
!     count interior cells, external=receive cells
!     count cell only once on external cells
!
      ip1=0
      DO ip=1,nnp
         ncell_fluid=ncell_fluid1(ip)
         ALLOCATE(jperm(ncell_fluid))
         cinter=0
         DO i=1,nelem
            flag(i)=0
            jperms(i)=0
         ENDDO
         DO jp=1,nnp
            irecv_cnt(jp)=0
         ENDDO
         iptr=0
         DO jj=ia(ip),ia(ip+1)-1 
            ie=ja(jj)
            DO j=ia_sub(ie),ia_sub(ie+1)-1
               k=ja_sub(j)
               jp=celem(k)
               IF(jp.eq.ip) cycle
!     cext cells
               DO j1=j,ia_sub(ie+1)-1
                  k=ja_sub(j1)
                  jp=celem(k)
                  IF(jp.eq.ip) cycle
                  IF(flag(k).eq.0) then
                     flag(k)=1
                     iptr=iptr+1
                     flagt(iptr)=k
                     irecv_cnt(jp)=irecv_cnt(jp)+1
                  ENDIF
               ENDDO
               GOTO 100
            ENDDO
!     cintr cells
!     get the mapping for interior cell
            cinter=cinter+1
            jperm(cinter)=ie
            jperms(ie)=cinter
100         CONTINUE 
            cinterp(ip)=cinter
         ENDDO
         cext0=0
         DO jp=1,nnp
            cext0=cext0+irecv_cnt(jp)
         ENDDO
!     minimize the zeroing of flag via flagt
         DO i=1,cext0
            k=flagt(i)
            flag(k)=0
         ENDDO
         IF(lev_typedt.eq.0) THEN
!...........No reordering
            DO i=1,cinter
               perm(ip1+i)=jperm(i)
               permi(jperm(i))=ip1+i
!              write(*,*) ip,ip1+i,perm(ip1+i)
            ENDDO
!           write(*,*) 
         ELSE
!...........reordering Metis,cuthill
!...........Build local a    
!           write(*,*) (jperm(jj),jj=1,cinter)
            ALLOCATE(ia2(cinter+1))
            ALLOCATE(perm2(cinter),permi2(cinter))
            ia2(1)=1
            maxmt2=0
            DO i=1,cinter
               ie=jperm(i)
!              write(*,*) (ja_sub(jj),jj=ia_sub(ie),ia_sub(ie+1)-1)
               do jj=ia_sub(ie),ia_sub(ie+1)-1
                  j=ja_sub(jj)
                  j1=jperms(j)
                  if(j1.eq.0) cycle
                  maxmt2=maxmt2+1
!                 write(*,*) '=>',j,j1
               enddo
               ia2(i+1)=maxmt2+1
            ENDDO
!              write(*,*) '*>',ip,maxmt2
            ALLOCATE(ja2(maxmt2))
            ip2=1
            DO i=1,cinter
               ie=jperm(i)
               do jj=ia_sub(ie),ia_sub(ie+1)-1
                  j=ja_sub(jj)
                  j1=jperms(j)
                  if(j1.eq.0) cycle
                  ja2(ip2)=j1
                  ip2=ip2+1
               enddo
            ENDDO
            DO i=1,cinter
!              write(*,*) i,(ja2(jj),jj=ia2(i),ia2(i+1)-1)
            ENDDO
            IF(lev_typedt.eq.1) THEN
               vwgt(:)=1
               CALL metis_setdefaultoptions(moptions)
               moptions(17)=1
               CALL METIS_NodeND(cinter,ia2,ja2,vwgt,moptions,perm2,permi2)
            ELSEIF(lev_typedt.eq.2) THEN
               CALL cuthill(cinter,ia2,ja2,perm2,permi2)
            ENDIF
            DO i=1,cinter
!              write(*,*) i,perm2(i),permi2(i)
            ENDDO
            DO i=1,cinter
               j1=perm2(i)
               perm(ip1+i)=jperm(j1)
               permi(jperm(j1))=ip1+i
!              write(*,*) i,perm(ip1+i),permi(jperm(j1))
            ENDDO
!              write(*,*) '======'
            DEALLOCATE(ia2,ja2)
            DEALLOCATE(perm2,permi2)
         ENDIF
         DEALLOCATE(jperm)
         ip1=ip1+cinter
      ENDDO
!     stop 10
!
!     count remote=send cells
!     count cell only once on remote cells for each ip
!
      ip1=0
      cinter=0
      DO ip=1,nnp
         cinter=cinter+cinterp(ip)
      enddo
         ip1=cinter
!     write(*,*) 'cinter',cinter
      DO jp=1,nnp
         ncell_fluid=ncell_fluid1(jp)
!        cinter=cinterp(jp)
         iptr=0
         DO ip=1,nnp
            jsend_cnt(ip)=0
            IF(jp.ne.ip) THEN
               DO jj=ia(ip),ia(ip+1)-1 
                  ie=ja(jj)
!..............Bypass scan of interior points
                  IF(jperms(ie).eq.0) THEN
                     DO j=ia_sub(ie),ia_sub(ie+1)-1
                        k=ja_sub(j)
                        j0=celem(k)
                        IF(j0.eq.jp) THEN
                           IF(flag(k).eq.0) then
                              flag(k)=1
                              jsend_cnt(ip)=jsend_cnt(ip)+1
                              flagt(iptr+jsend_cnt(ip))=k
                           ENDIF
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO
!...........Minimize the zeroing of flag via flagt
               DO i=1,jsend_cnt(ip)
                  k=flagt(iptr+i)
                  flag(k)=0
               ENDDO
               iptr=iptr+jsend_cnt(ip)
            ENDIF
         ENDDO
!
!     send cells
!
!     remove zero count
         niut=0
         DO ip=1,nnp
            IF(jsend_cnt(ip).ne.0) niut=niut+1
         ENDDO
         ALLOCATE(iut(niut),si(niut+1))
         si(1)=1
         j=0
         DO ip=1,nnp
            IF(jsend_cnt(ip).ne.0) THEN
               j=j+1
!              iut(j)=jp
!              iutjp(jp)=j
               si(j+1)=jsend_cnt(ip)
            ENDIF
         ENDDO
         DO j=1,niut
            si(j+1)=si(j+1)+si(j)
         ENDDO
!
         iut=iut-1
!
!     get all neighbors counted once jsend
!
         ALLOCATE(jsend(si(niut+1)-1))
         DO i=1,si(niut+1)-1
            jsend(i)=flagt(i)
        ENDDO
!
!     get the mapping for remote cells counted once.
!
         DO i=1,si(niut+1)-1
            k1=jsend(i)
            IF(flag(k1).eq.0) then
               k=k+1
               ip1=ip1+1
               flag(k1)=1
               perm(ip1)=k1
               permi(k1)=ip1
!              write(*,*) jp,ip1,k1
            ENDIF
         ENDDO
!        ip1=ip1+ncell_fluid
         DEALLOCATE(iut,si)
         DEALLOCATE(jsend)
      ENDDO
!     do i=1,nelem
!        write(*,*) 'perm',i,perm(i),permi(i)
!     enddo
!      stop 10
!
!DEC$ELSE
      IF(myrank.eq.0)WRITE(*,*)'domain_decomposition is deactivated due to false metis_flag!!!'
      PAUSE
      STOP
!DEC$ENDIF
      END SUBROUTINE reorder_domain
!
      SUBROUTINE final_r(n,maxmt1,isize,                 &
                         ia,iend,ja,                     &
                         ia1,ju,ja1,                     &
                         perm,permi,indexi)
!                        perm,permi,indexi,nnp,cinterp)
!
      IMPLICIT NONE
!              
!.....Input
      INTEGER :: n,maxmt1,isize,nnp
      INTEGER :: ia(n+1),iend(n),ja(maxmt1)
      INTEGER :: perm(n),permi(n)
!     INTEGER :: cinterp(nnp)
!.....Output
      INTEGER :: ia1(n+1),ju(n),ja1(isize)
      INTEGER :: indexi(isize)
!.....local variables
      integer index(n)
!
      integer i,i1,j,j1,jj,ip,nn
!
!        write(*,*) 'i==>',maxmt1,isize
!        stop 99
      ia1(1)=1
      ip=1
      do i=1,n
         i1=perm(i)
!        write(*,*) 'i==>',i,i1
!        write(*,400) (ja(jj),jj=ia(i1),iend(i1))
         do jj=ia(i1),iend(i1)
            j=ja(jj)
            j1=permi(j)
            ja1(ip)=j1
!           index(jj-ia(i1)+1)=jj-ia(i1)+1
            ip=ip+1
         enddo
         ia1(i+1)=ip
!        write(*,400) (ja1(jj),jj=ia1(i),ia1(i+1)-1)
400      format(20(i7,2x))
!
         jj=ia1(i)
         nn=ia1(i+1)-ia1(i)
       do j=1,nn
          index(j)=j
       enddo
         call sortx_i(ja1(jj),index,nn)
         DO j=1,nn
            j1=index(j)
            indexi(jj+j1-1)=j
         ENDDO
!        write(*,400) (ja1(jj),jj=ia1(i),ia1(i+1)-1)
!        stop 10
      enddo
      do i=1,n
         do j=ia1(i),ia1(i+1)-1
            IF(ja1(j).eq.i) THEN
               ju(i)=j
               exit
            ENDIF
         enddo
      enddo
!     cintertot=0
!     do i=1,nnp
!        cintertot=cintertot+cinterp(i) 
!     enddo
!     tot1=0
!     do i=cintertot+1,n
!        do jj=ia1(i),ia1(i+1)-1
!           j=ja1(jj)
!           if(j.gt.cintertot) tot1=tot1+1
!        enddo
!     enddo
!        write(*,400) ia1(n+1)-1,tot1
!        write(*,*) 100.d0*tot1/(ia1(n+1)-1)
!     do i=1,n
!        write(*,400) i,(ja1(jj),jj=ia1(i),ia1(i+1)-1)
!     enddo
!     stop 98
!
      END SUBROUTINE final_r
!
      SUBROUTINE copya1_r(n,                  &
                         au,ia_a,ja_a,       &
                         ia_r,ja_r,a_r, &
                         perm,permi,index)
!
!
      IMPLICIT NONE
!
!.....Input
      INTEGER :: n
      INTEGER :: ia_a(n+1),ja_a(*)
      INTEGER :: ia_r(n+1),ja_r(*)
      INTEGER :: perm(n),permi(n),index(*)
      REAL(8) :: au(*)
!.....Output
      REAL(8) :: a_r(*)
!.....local variables
      INTEGER :: i,j,jj,i1,j1,ip,k
      INTEGER :: l,j2

!.....Local arrays
      INTEGER :: ja2(n)
      REAL(8) :: a1(n)
!
      ip=1
      DO i=1,n
         i1=perm(i)
         l=ia_a(i1+1)-ia_a(i1)
         DO jj=ia_a(i1),ia_a(i1+1)-1
            j=ja_a(jj)
            j1=permi(j)
            k=index(ip)
            a1(k)=au(jj)
            ja2(k)=j1
            ip=ip+1
         ENDDO
         j1=1
         DO jj=ia_r(i),ia_r(i+1)-1
            IF(j1.gt.l) THEN
               DO j2=jj,ia_r(i+1)-1
                  a_r(j2)=0.d0
               ENDDO
               exit
            ENDIF
            IF(ja_r(jj).eq.ja2(j1)) THEN
               a_r(jj)=a1(j1)
               j1=j1+1
            ELSE
               a_r(jj)=0.d0
            ENDIF
         ENDDO
      ENDDO
!
      END SUBROUTINE copya1_r
!
      subroutine dgefs2s(n,maxmt,maxmt1, &
                         ia,ja,          &
                         ia1,ja1,idiags)
!
      implicit none
!     input
      integer n,maxmt,maxmt1
      integer ia(n+1),ja(maxmt)
!     output
      integer ia1(n+1),ja1(maxmt1)
      integer idiags(n)
!     local variables
      integer i,j,jf,jj,kk
      integer i0,i1,ip1,ip2,ip,k1,k2
!     local arays
      integer jat(n,0:1)
!
      ia1(1)=1
      do i=1,n
!           write(*,*) '++',i
! copy row i of a1 to buff at(,0)  get jf lenght of ja
         i0=0
         i1=1
!        idiags(i)=ia1(i)
         jf=0
         do jj=ia(i),ia(i+1)-1
            jf=jf+1
            jat(jf,i0)=ja(jj)
         enddo
!        write(*,*) (ja(jj),jj=ia(i),ia(i+1)-1)
         jj=1
100      continue
         if(jj.gt.jf) goto 110
! get t1
         j=jat(jj,i0)
         if(j.lt.i) then
            jat(jj,i1)=j
            jat(jj,i0)=j
         elseif(j.eq.i) then
! treat pivot            
            idiags(i)=ia1(i)+jj-1
            jat(jj,i0)=j
            goto 110
         endif
! compute one row only
            ip1=idiags(j)+1
            ip2=jj+1
            ip=jj
200         continue
            if(ip1.gt.ia1(j+1)-1) then
               do kk=ip2,jf
                  ip=ip+1
                  jat(ip,i1)=jat(kk,i0)
               enddo
               goto 210
            endif
            if(ip2.gt.jf) then
               do kk=ip1,ia1(j+1)-1
                 ip=ip+1
                 jat(ip,i1)=ja1(kk)
               enddo
               goto 210
            endif
            k1=ja1(ip1)
            k2=jat(ip2,i0)
            if(k1.lt.k2) then
               ip=ip+1
               jat(ip,i1)=k1
               ip1=ip1+1
            elseif(k1.gt.k2) then
               ip=ip+1
               jat(ip,i1)=k2
               ip2=ip2+1
            else
               ip=ip+1
               jat(ip,i1)=k2
               ip1=ip1+1
               ip2=ip2+1
            endif
            goto 200
210         continue
            i0=mod(i0+1,2)
            i1=mod(i1+1,2)
            jj=jj+1
            jf=ip
         goto 100
110      continue
!        write(*,*) i,ia1(i),jf,ia1(i)+jf-1
         i1=ia1(i)
         if(i1+jf-1.gt.maxmt1) then
!           write(*,*) 'error alloc jas',i,maxmt1,ia1(i+1)-1
            write(*,*) 'error alloc jas',i1,jf,i1+jf-1,maxmt1
            stop 98
         endif
         do jj=1,jf
            ja1(i1)=jat(jj,i0)
            i1=i1+1
         enddo
         ia1(i+1)=i1
      enddo
!
      END SUBROUTINE dgefs2s
      subroutine dgefs2s_lev(n,maxmt,maxmt1,          &
                             ia,ja,                   &
                             ia1,ja1,idiags,lev,lev0)
!
      implicit none
!     input
      integer n,lev0
      integer maxmt,maxmt1
      integer ia(n+1),ja(maxmt)
!     output
      integer ia1(n+1),ja1(maxmt1),lev(maxmt1)
      integer idiags(n)
!     local variables
      integer i,j,jf,jj,kk,levik
      integer i0,i1,ip1,ip2,ip,k1,k2
!     local arays
      integer jat(n,0:1)
      integer levt(n,0:1)
!
      ia1(1)=1
      do i=1,n
! copy row i of a1 to buff at(,0)  get jf lenght of ja
         i0=0
         i1=1
!        idiags(i)=ia1(i)
         jf=0
         do jj=ia(i),ia(i+1)-1
            jf=jf+1
            jat(jf,i0)=ja(jj)
            levt(jf,i0)=0
         enddo
!         
         jj=1
100      continue
         if(jj.gt.jf) goto 110
! get t1
         j=jat(jj,i0)
         if(j.lt.i) then
            jat(jj,i1)=j
            jat(jj,i0)=j
            levik=levt(jj,i0)
            levt(jj,i1)=levik
         elseif(j.eq.i) then
! treat pivot            
            idiags(i)=ia1(i)+jj-1
            jat(jj,i0)=j
            goto 110
         endif
! compute one row only
            ip1=idiags(j)+1
            ip2=jj+1
            ip=jj
200         continue
            if(ip1.gt.ia1(j+1)-1) then
               do kk=ip2,jf
                  ip=ip+1
                  jat(ip,i1)=jat(kk,i0)
                  levt(ip,i1)=levt(kk,i0)
               enddo
               goto 210
            endif
            if(ip2.gt.jf) then
               do kk=ip1,ia1(j+1)-1
                 if(levik+lev(kk)+1.le.lev0) then
                 ip=ip+1
                 jat(ip,i1)=ja1(kk)
                 levt(ip,i1)=levik+lev(kk)+1
                 endif
               enddo
               goto 210
            endif
            k1=ja1(ip1)
            k2=jat(ip2,i0)
            if(k1.lt.k2) then
               if(levik+lev(ip1)+1.le.lev0) then
               ip=ip+1
               jat(ip,i1)=k1
               levt(ip,i1)=levik+lev(ip1)+1
               endif
               ip1=ip1+1
            elseif(k1.gt.k2) then
               ip=ip+1
               jat(ip,i1)=k2
               levt(ip,i1)=levt(ip2,i0)
               ip2=ip2+1
            else
               if(min(levt(ip2,i0),levik+lev(ip1)+1).le.lev0) then
               ip=ip+1
               jat(ip,i1)=k2
               levt(ip,i1)=min(levt(ip2,i0),levik+lev(ip1)+1)
               endif
               ip1=ip1+1
               ip2=ip2+1
            endif
            goto 200
210         continue
            i0=mod(i0+1,2)
            i1=mod(i1+1,2)
            jj=jj+1
            jf=ip
         goto 100
110      continue
         i1=ia1(i)
         if(i1+jf-1.gt.maxmt1) then
            write(*,*) 'error alloc0jas',i,maxmt1,ia1(i+1)-1
            stop 99
         endif
         do jj=1,jf
            ja1(i1)=jat(jj,i0)
            lev(i1)=levt(jj,i0)
            i1=i1+1
         enddo
         ia1(i+1)=i1
      enddo
!
      END SUBROUTINE dgefs2s_lev
      subroutine ddgef2(a,lda,n)
      implicit none
      integer lda,n
      real*8 a(lda,n)
      integer i,j,k
      real*8 pivotx
!
      do i=1,n
!     if(mod(i,100).eq.0) write(*,*) i
         do k=1,i-1
           if(a(i,k).ne.0.d0) then
            a(i,k)=a(i,k)*a(k,k)
            do j=k+1,i-1
               a(i,j)=a(i,j)-a(i,k)*a(k,j)
            enddo
           endif
         enddo

            j=i
         do k=1,i-1
!           a(i,i)=a(i,i)-a(i,k)*a(k,i)
         enddo
         do k=1,i-1
           if(a(i,k).ne.0.d0) then
            do j=i,n
               a(i,j)=a(i,j)-a(i,k)*a(k,j)
            enddo
           endif
         enddo
         pivotx=1.d0/a(i,i)
         a(i,i)=pivotx
      enddo
!
      return
      end
       subroutine ddges(a,lda,n,b)
       implicit none
!
       integer lda,n
       real*8  a(lda,*)
       real*8  b(*)
!
       integer i,j
       real*8  temp

!------
          do i=2,n
             temp=b(i)
             do j=1,i-1
                temp=temp-a(i,j)*b(j)
             enddo
             b(i)=temp
          enddo

!         do i=1,n
!            i=n
!            write(6,*) i,b(i)
!         enddo
!------
          do i=n,1,-1
             temp=b(i)
             do j=i+1,n
                temp=temp-a(i,j)*b(j)
             enddo
             b(i)=temp*a(i,i)
          enddo

!         write(6,*) 'after'
!         do i=1,n
!            i=n
!            write(6,*) i,b(i)
!         enddo
!
       return
       end
