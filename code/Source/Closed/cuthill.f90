      SUBROUTINE cuthill(n,ia,ja,perm,permi)
!
      IMPLICIT NONE
!
!.....Input
      INTEGER :: n
      INTEGER :: ia(*),ja(*)
!.....Output
      INTEGER :: perm(*),permi(*)
!.....Local variables
      INTEGER :: i,i1
!
      CALL genrcm(n,perm,ia,ja)
      DO i=1,n
         i1=perm(i)
         permi(i1)=i
      ENDDO
!
      END SUBROUTINE cuthill
!
      SUBROUTINE genrcm(n,perm,ia,ja)
!
      IMPLICIT NONE
!
!.....Input
      INTEGER :: n
      INTEGER :: ia(*),ja(*)
!.....Output
      INTEGER :: perm(*)
!.....Local variables
      INTEGER :: i,num,nlvl
      INTEGER :: root,nums
!.....Local arrays
      INTEGER :: mask(n)
      INTEGER :: ias(n+1)
!
      DO i=1,n
         mask(i)=1
      ENDDO
      num=1
      DO i=1,n
         IF(mask(i).ne.0) then
            root=i
            CALL fnroot(root,n,nlvl,ias,perm(num),ia,ja)
            CALL rcm(root,n,nlvl,mask,nums,ias,perm(num),ia,ja) 
            num=num+nums
            IF(num.gt.n) RETURN
         ENDIF
      ENDDO
!
      END SUBROUTINE genrcm
!
      SUBROUTINE fnroot(root,n,nlvl,ias,jas,ia,ja)
!
      IMPLICIT NONE
!
!.....Input
      INTEGER :: n,nlvl
      INTEGER :: ias(*),jas(*)
      INTEGER :: ia(*),ja(*)
!.....Output
      INTEGER :: root
!.....Local variables
      INTEGER :: nums
      INTEGER :: j,jj,js,k,mindeg,node,ndeg
      INTEGER :: numlvl
!
      CALL rootls(root,n,nlvl,ias,jas,ia,ja)
      nums=ias(nlvl+1)-1
! 
      IF(nlvl.eq.1 .or. nlvl.eq.nums) RETURN
!
100   CONTINUE
      js=ias(nlvl)
      mindeg=nums
      root=jas(js)
      DO jj=ias(nlvl),ias(nlvl+1)-1
         node=jas(jj)
         ndeg=0
         DO k=ia(node),ia(node+1)-1
            j=ja(k)
            IF(j.ne.node) ndeg=ndeg+1
         ENDDO
         IF(ndeg.lt.mindeg) then
            root=node
            mindeg=ndeg
         ENDIF
      ENDDO
      CALL rootls(root,n,numlvl,ias,jas,ia,ja)
      IF(numlvl.le.nlvl) goto 200
        nlvl=numlvl
        IF(nlvl.lt.nums) goto 100
200   CONTINUE
!
      END SUBROUTINE fnroot
!
      SUBROUTINE rcm(root,n,nlvl,mask,nums,ias,perm,ia,ja)
!
      IMPLICIT NONE
!
!.....Input
      INTEGER :: n,nlvl,root
      INTEGER :: ia(*),ja(*)
!.....Output
      INTEGER :: mask(n)
      INTEGER :: perm(*)
      INTEGER :: ias(*)
!.....Local variables
!.....Local arrays
      INTEGER :: ns
      INTEGER :: fnbr
      INTEGER :: nums,deg(n)
      INTEGER :: i,jj,node,nbr,l,lperm
      INTEGER :: lbegin,lvlend,lnbr
      INTEGER,DIMENSION(:),ALLOCATABLE :: sort
!
      nums=ias(nlvl+1)-1
      IF(nums.lt.1) return
      CALL degree(nlvl,deg,ias,perm,ia,ja,ns)      
      ALLOCATE(sort(ns))
!
      mask(root)=0
      lvlend=0
      lnbr=1
100   CONTINUE
      lbegin=lvlend+1
      lvlend=lnbr
      DO i=lbegin,lvlend
         node=perm(i)
         fnbr=lnbr+1
         DO jj=ia(node),ia(node+1)-1
            nbr=ja(jj)
            IF(nbr.eq.node) cycle
            IF(mask(nbr).ne.0) then
               lnbr=lnbr+1
               mask(nbr)=0
               perm(lnbr)=nbr
               sort(lnbr-fnbr+1)=deg(nbr)
            ENDIF
         ENDDO
         CALL sortx_i(sort,perm(fnbr),lnbr-fnbr+1)
      ENDDO
      IF(lnbr.gt.lvlend) goto 100
!
      DEALLOCATE(sort)
!reverse
      l=nums
      DO i=1,nums/2
         lperm=perm(l)
         perm(l)=perm(i)
         perm(i)=lperm
         l=l-1
      ENDDO
!
      END SUBROUTINE rcm
! 
      SUBROUTINE rootls(root,n,nlvl,ias,jas,ia,ja)
!
      IMPLICIT NONE
!
!.....Input
      INTEGER :: n,root
      INTEGER :: ia(*),ja(*)
!.....Output
      INTEGER :: nlvl
      INTEGER :: ias(*),jas(*)
!.....Local variables
      INTEGER :: i,j,jj,i1,node
!.....Local arrays
      INTEGER :: mask(n)
!
      DO i=1,n
         mask(i)=1
      ENDDO
!
      ias(1)=1
      jas(1)=root
      mask(root)=0
      nlvl=0
      i1=2
200   CONTINUE
      nlvl=nlvl+1
      ias(nlvl+1)=i1
      DO i=ias(nlvl),ias(nlvl+1)-1
         node=jas(i)
         DO jj=ia(node),ia(node+1)-1
            j=ja(jj)
            IF(j.eq.node) cycle
            IF(mask(j).ne.0) THEN
               jas(i1)=j
               mask(j)=0
               i1=i1+1
            ENDIF
         ENDDO
      ENDDO
      IF(i1.gt.ias(nlvl+1)) goto 200
!
      END SUBROUTINE rootls
!
      SUBROUTINE degree(nlvl,deg,ias,jas,ia,ja,ns)
!
      IMPLICIT NONE
!
!.....Input
      INTEGER :: nlvl
      INTEGER :: ias(*),jas(*)
      INTEGER :: ia(*),ja(*)
!.....Output
      INTEGER :: ns
      INTEGER :: deg(*)
!.....Local variables
      INTEGER :: i,j,k,jj,l,node,ideg
      INTEGER :: nums
!
      ns=0
      DO i=1,nlvl
         DO jj=ias(i),ias(i+1)-1
            node=jas(jj)
            l=ia(node+1)-ia(node)
            ns=max(ns,l)
            ideg=0
            DO k=ia(node),ia(node+1)-1
               j=ja(k)
               IF(j.ne.node) ideg=ideg+1
            ENDDO
            deg(node)=ideg
         ENDDO
      ENDDO
      nums=ias(nlvl+1)-1
!
      END SUBROUTINE degree
