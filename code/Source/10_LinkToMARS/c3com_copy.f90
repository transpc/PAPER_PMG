!
      SUBROUTINE c3com_copy_C2S
!
!     This routine copies c3com data
!
      USE Vol_data
      USE c3com_cupid
      IMPLICIT none
!      
      INCLUDE 'c3com.h' 
      INCLUDE 'c3com_space.h' 
!
!====================      
!.....INTEGER parameter
!====================      
!!!      j1cupid(:) = i3cupid(:)
!
!====================      
!.....REAL parameter
        !====================      
        s3vpp(:)=cell%pps(i3cupid(:))
        s3ngpp(:)=cell%p(i3cupid(:))-cell%pps(i3cupid(:))
        s3vt(:)=cell%tg(i3cupid(:))
        s3lt(:)=cell%tl(i3cupid(:))
        s3dt(:)=cell%td(i3cupid(:))
        s3ngmf(:,:)=0.0d0
      
        s3pa(:)=c3pa(1,:)
        s3ug(:)=c3ug(1,:)
        s3uf(:)=c3uf(1,:)
        s3al(:)=c3al(1,:) 
        s3rhof(:)=c3rhof(1,:)
        s3rhog(:)=c3rhog(1,:)
        s3brn(:)=c3brn(1,:)
        s3arxq(:)=c3arxq(1,:)
        s3vpgno(:)=0.0d0   !MCC-jjj-NEXT
      
!      s3ngmf(:,1)=qn_cell(i3cupid(:),3)
!      s3ngmf(:,2)=qn_cell(i3cupid(:),1)
!      s3ngmf(:,3)=qn_cell(i3cupid(:),2)
!      s3ngmf(:,4)=qn_cell(i3cupid(:),4)
!      s3ngmf(:,5)=qn_cell(i3cupid(:),5)
!      s3ngmf(:,6)=qn_cell(i3cupid(:),7)
!      s3ngmf(:,7)=qn_cell(i3cupid(:),6)
!      s3ngmf(:,8)=qn_cell(i3cupid(:),8)
!      s3ngmf(:,9)=0.0d0
!      s3ngmf(:,10)=0.0d0
      
      ! SPACE         CUPID
      ! 0 : N2        1 : He
      ! 1 : He        2 : H2
      ! 2 : H2        3 : N2
      ! 3 : Kr        4 : Kr
      ! 4 : Xe        5 : Xe
      ! 5 : Ar        6 : Air
      ! 6 : Air       7 : Ar
      ! 7 : SF6       8 : SF6
      ! 8 : user1
      ! 9 : user2

      RETURN 
      END SUBROUTINE 

!======================================================================

!
      SUBROUTINE c3com_copy_S2C
!
!     This routine copies c3com data
!
      USE c3com_cupid

      IMPLICIT none
!      
      INCLUDE 'c3com.h' 
      INCLUDE 'c3com_space.h' 
!
!
!.....LOGICAL parameters
!
      IF(jflag_relap.eq.1)then
         flag_relap=.true.
      ELSE
         flag_relap=.false.
      ENDIF
      
!!!      IF(jflag_cobra.eq.1)then
!!!         flag_cobra=.true.
!!!      ELSE
!!!         flag_cobra=.false.
!!!      ENDIF
      
      IF(jflag_stop.eq.1)then
         flag_stop=.true.
      ELSE
         flag_stop=.false.
      ENDIF
      IF(jflag_bd.eq.1)then
         flag_bd=.true.
      ELSE
         flag_bd=.false.
      ENDIF
      IF(jflag_wp.eq.1)then
         flag_wp=.true.
      ELSE
         flag_wp=.false.
      ENDIF
      IF(jflag_bd2.eq.1)then
         flag_bd2=.true.
      ELSE
         flag_bd2=.false.
      ENDIF
      IF(s3sdbvolpack.eq.1)then
         sdbvolpack=.true.
      ELSE
         sdbvolpack=.false.
      ENDIF
!      IF(s3sdbvolpacko(:).eq.1)then
!         sdbvolpacko(1,:)=.true.
!      ELSE
!         sdbvolpacko(1,:)=.false.
!      ENDIF
!      IF(s3c3pack(:,:).eq.1)then
!         c3pack(:,:)=.true.
!      ELSE
!         c3pack(:,:)=.false.
!      ENDIF
!      IF(s3overcorrection(:).eq.1)then
!         overcorrection(1,:)=.true.
!      ELSE
!         overcorrection(1,:)=.false.
!      ENDIF
!====================      
!.....REAL parameter
!====================
!!!         s3pa(:)=c3pa(1,:)
!!!         s3ug(:)=c3ug(1,:)
!!!         s3uf(:)=c3uf(1,:)
!!!         s3al(:)=c3al(1,:) 
!!!         s3rhof(:)=c3rhof(1,:)
!!!         s3rhog(:)=c3rhog(1,:)
!!!         s3brn(:)=c3brn(1,:)
!!!         s3arxq(:)=c3arxq(1,:)
!!!         s3vpgno(:)=0.0d0   !MCC-jjj-NEXT
            
      dt_super           = s3dt_super    
      dt_relap           = s3dt_relap    
!!!      dt_cobra           = s3dt_cobra    
      c3vg(1,:)         = s3vg(:)      
      c3vl(1,:)         = s3vl(:)      
      c3alphg(1,:)      = s3alphg(:)   
      c3alphf(1,:)      = s3alphf(:)   
      c3betag(1,:)      = s3betag(:)   
      c3betaf(1,:)      = s3betaf(:)   
      c3xi(1,:)         = s3xi(:)      
      c3yeta(1,:,:)     = s3yeta(:,:)
      c3area(1,:)       = s3area(:)    
      c3rtp(1,:,:)      = s3rtp(:,:)  
      c3delp(1,:)       = s3delp(:)    
      c3dpmt(1,:,:)     = s3dpmt(:,:)
!!!      c3pa(1,:)         = s3pa(:)    
!!!      c3uf(1,:)         = s3uf(:)    
!!!      c3ug(1,:)         = s3ug(:)    
!!!      c3al(1,:)         = s3al(:)    
!!!      c3arxq(1,:)       = s3arxq(:)  
!      c3vvl(1,:)        =      
!      c3vvg(1,:)        =      
      c3odr(:)           = s3odr(:)
      c3dpm(:)           = s3dpm(:)
!!!      c3rhof(1,:)       = s3rhof(:) 
!!!      c3rhog(1,:)       = s3rhog(:) 
      r3ent(1,:)        = s3ent(:)  
      r3liq(1,:)        = s3liq(:)  
      cc3rtp(1,:,:)    = ss3rtp(:,:) 
!!!      c3brn(1,:)        = s3brn(:)
!!!      c3vpgno(1,:)      = s3vpgno(:)
!      packdvf(1,:)     = 
!      aloold(100,100)   =  s3aloold(100,10)
      packfactor         = s3packfactor
!                        
!====================      
!.....INTEGER parameter       
!====================      
!      ibeg50            = 
      i_where            = j_where
      i3bcn(1,:)        = j3bcn(:)
!      n_repet           = 
!      nstep_c           = 
!      nstep_r           = 
!      iend50            = 
      i3nodr             = j3nodr
      i3nic(:)           = j3nic(:)
      i3chan(1,:)       = j3chan(:)
      i3cell(1,:)       = j3cell(:)
      i3mode(1,:)       = j3mode(:)
      i3line(:)         = j3line(:)
      i3modet(1,:,:)    = j3modet(:,:)
      i3dir(1,:)        = j3dir(:)  
      i1Cvoln(1,:)      = j1Cvoln(:)
      i1Cvodn(1,:)      = j1Cvodn(:)
      i1Cvndx(1,:)      = j1Cvndx(:)
      i1jndx(1,:)       = j1jndx(:) 
      i1Rvodn(1,:)      = j1Rvodn(:)
      i1Rvndx(1,:)      = j1Rvndx(:)
      hpackflag(:,:) = s3hpackflag(:,:)
      vpackflag(:,:) = s3vpackflag(:,:)
      i3cupid(:)    = j1cupid(:)
!
        
      RETURN 
      END SUBROUTINE 
