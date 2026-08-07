! 
!.....................................................................  
!     common variables for the cobra/relap5 code integration             
!                                                                        
!     written by S.Y.Lee, 1991 on Apollo DN10000.                        
!     modified by J.J.Jeong, 08/07/1996 for HP 735                       
!                                                                       
!     restructured by J.J.Jeong, April 1997 for HP & PC                  
!     array size changed by J.J.Jeong, May 1998   (5,6) --> (1,20)
!     array size changed by J.J.Jeong, March 1999 (1,20) --> (1,40)
!     array size changed by I.K.Park, January 2012 (1,40) --> (1,72) !mcc-pik
!     variable related to 3_D packing is added by h.g.lim june 2002
!    : sdbvolpack, sdbvolpacko,packdvf
!
!     Note: indxsp(10,72) & i1nic(10,72) in c3connct.h were also changed.
!.....................................................................  
! 
!

      COMMON/c3com_space/        &
      s3dt_super               , & !final time step size
      s3dt_relap               , & !1-D code time step size
      s3dt_cobra               , & !3-D code time step size
      s3vg                     , & !explicit vapour velocity at face between 1-D and 3-D regions
      s3vl                     , & !  
      s3alphg                  , & !explicit vapour velocity = c3vg
      s3alphf                  , & !  
      s3betag                  , & !coefficient of pressure variation in 1-D cell
      s3betaf                  , & !  
      s3xi                     , & !pressure variation in 1-D cell
      s3yeta                   , & !coefficient of pressure variation in 3-D cell
      s3area                   , & !face area at interface between 1-D and 3-D regions
      s3rtp                    , & !interface data(void fractions, NCG quality, densities, energies, P...) donor properties
      s3delp                   , & !pressure variation in 3-D cell (from 3-D pressure eq.)
      s3brn                    , & !boron concentration in 3-D cell (from 3-D boron eq.)
      s3vpgno                  , & !vapor generation in 3-D cell (from 3-D calculation result)
      s3dpmt                   , & !used in 3-D
      s3pa                     , & !pressure in 3-D cell for 3D cell TH properties (from 3-D boron eq.)
      s3uf                     , & !internal energy of liquid in 3-D cell TH properties
      s3ug                     , & !  
      s3al                     , & !void fraction in 3-D cell
      s3arxq                   , & !NCG quality in 3-D cell
      s3odr                    , & !variable of 3-D cell eq.
      s3dpm                    , & !variable of 3-D cell eq.
      s3aloold                 , & !
      s3packfactor             , & !      
      s3rhof                   , & !liquid density at 3-D 
      s3rhog                   , & !  
      s3ent                    , & !liquid partition
      s3liq                    , & !liquid partition
      ss3rtp                   , & !old time c3rtp
	                               !!!!!!!!!!!!!!!!!!!!!!!!
	                               !!!!! ADDED by C.W.Choi
	                               !!!!!!!!!!!!!!!!!!!!!!!!
	  s3vpp                    , & !Vapor partial pressure in 3D cell (BC in SPACE)
	  s3ngpp                   , & !Noncondensible partial pressure in 3D cell (BC in SPACE)
	  s3vt                     , & !Vapor temperature in 3D cell (BC in SPACE)
	  s3lt                     , & !Liquid temperature in 3D cell (BC in SPACE)
	  s3dt                     , & !Droplet temperature in 3D cell (BC in SPACE)
	  s3ngmf                   , & !Noncondensible gas mass fractio in 3D cell (BC in SPACE)
	                               !!!!!!!!!!!!!!!!!!!!!!!!
	                               !!!!! ADDED by C.W.Choi
	                               !!!!!!!!!!!!!!!!!!!!!!!!
	  j1cupid                  , & !CUPID component number connected to SPACE
      j_where                  , & !step identifier being processed in transient routine
      jflag_relap              , & !indicator for calculation repetition for 1-D code
      jflag_cobra              , & !indicator for calculation repetition for 3-D code
      jflag_stop               , & !indicator for calculation stop
      jflag_bd                 , & !indicator for bad donor
      jflag_wp                 , & !indicator for water packing
      jflag_bd2                , & !indicator for bad donor after water packing correction
      j3bcn                    , & !3-D cell index, used in 3-D 
      j1Cvoln                  , & !3-D cell number at 1-D
      j3chan                   , & !     
      j3cell                   , & !     
      j3mode                   , & !     
      j3modet                  , & !  
      j3nodr                   , & !number of 1-D regions
      j3nic                    , & !number of connections in a 1-D region
      j3line                   , & !     
      j3dir                    , & !connection direction used in only 1-D, + from 3-D to 1-D
      j1Cvodn                  , & !3-D cell index at 1-D
      j1Cvndx                  , & !3-D cell index at 1-D = i1Cvodn
      j1jndx                   , & !connection face index in 1-D face block
      j1Rvodn                  , & !interface 1-D cell index in cell block = i1Rvndx
      j1Rvndx                  , & !     
      jndxsp                   , & ! 
      j1nic                    , & ! 
      s3sdbvolpack             , & !         
      s3sdbvolpacko            , & !     
      s3c3pack                 , & !
      s3vpackflag              , & !
      s3hpackflag              , & !
      s3overcorrection         , & !     
      j1max                    

      REAL(8) s3dt_super           ! final time step size
      REAL(8) s3dt_relap           ! 1-D code time step size
      REAL(8) s3dt_cobra           ! 3-D code time step size
      REAL(8) s3vg(72)             ! explicit vapour velocity at face between 1-D and 3-D regions
      REAL(8) s3vl(72)             ! 
      REAL(8) s3alphg(72)          ! explicit vapour velocity = c3vg
      REAL(8) s3alphf(72)          ! 
      REAL(8) s3betag(72)          ! coefficient of pressure variation in 1-D cell
      REAL(8) s3betaf(72)          ! 
      REAL(8) s3xi(72)             ! pressure variation in 1-D cell
      REAL(8) s3yeta(72,72)       ! coefficient of pressure variation in 3-D cell
      REAL(8) s3area(72)           ! face area at interface between 1-D and 3-D regions
      REAL(8) s3rtp(72,15)         ! interface data(void fractions, NCG quality, densities, energies, P...) donor properties
      REAL(8) s3delp(72)           ! pressure variation in 3-D cell (from 3-D pressure eq.)
      REAL(8) s3brn(72)            ! boron concentration in 3-D cell (from 3-D boron eq.)
      REAL(8) s3vpgno(72)          ! vapor generation in 3-D cell (from 3-D calculation result)
      REAL(8) s3dpmt(72,6)         ! used in 3-D
      REAL(8) s3pa(72)             ! pressure in 3-D cell for 3D cell TH properties (from 3-D boron eq.)
      REAL(8) s3uf(72)             ! internal energy of liquid in 3-D cell TH properties
      REAL(8) s3ug(72)             ! 
      REAL(8) s3al(72)             ! void fraction in 3-D cell
      REAL(8) s3arxq(72)           ! NCG quality in 3-D cell
      REAL(8) s3odr(6)             ! variable of 3-D cell eq.
      REAL(8) s3dpm(6)             ! variable of 3-D cell eq.
      REAL(8) s3aloold(100,100)     ! 
      REAL(8) s3packfactor         ! 
      REAL(8) s3rhof(72)           ! liquid density at 3-D 
      REAL(8) s3rhog(72)           ! 
      REAL(8) s3ent(72)            ! liquid partition
      REAL(8) s3liq(72)            ! liquid partition
      REAL(8) ss3rtp(72,15)        ! old time c3rtp
	                               !!!!!!!!!!!!!!!!!!!!!!!!
	                               !!!!! ADDED by C.W.Choi
	                               !!!!!!!!!!!!!!!!!!!!!!!!
	  REAL(8) s3vpp(72)            !Vapor partial pressure in 3D cell (BC in SPACE)
	  REAL(8) s3ngpp(72)           !Noncondensible partial pressure in 3D cell (BC in SPACE)
	  REAL(8) s3vt(72)             !Vapor temperature in 3D cell (BC in SPACE)
	  REAL(8) s3lt(72)             !Liquid temperature in 3D cell (BC in SPACE)
	  REAL(8) s3dt(72)             !Droplet temperature in 3D cell (BC in SPACE)
	  REAL(8) s3ngmf(72,10)        !Noncondensible gas mass fractio in 3D cell (BC in SPACE)
	                               !!!!!!!!!!!!!!!!!!!!!!!!
	                               !!!!! ADDED by C.W.Choi
	                               !!!!!!!!!!!!!!!!!!!!!!!!
	  INTEGER j1cupid(72)          !CUPID component number connected to SPACE
      INTEGER j_where              ! step identifier being processed in transient routine
      INTEGER jflag_relap          ! indicator for calculation repetition for 1-D code
      INTEGER jflag_cobra          ! indicator for calculation repetition for 3-D code
      INTEGER jflag_stop           ! indicator for calculation stop
      INTEGER jflag_bd             ! indicator for bad donor
      INTEGER jflag_wp             ! indicator for water packing
      INTEGER jflag_bd2            ! indicator for bad donor after water packing correction
      INTEGER j3bcn(72)            ! 3-D cell index, used in 3-D 
      INTEGER j1Cvoln(72)          ! 3-D cell number at 1-D
      INTEGER j3chan(72)           !
      INTEGER j3cell(72)           !
      INTEGER j3mode(72)           !
      INTEGER j3modet(72,3)        !
      INTEGER j3nodr               ! number of 1-D regions
      INTEGER j3nic(3)             ! number of connections in a 1-D region
      INTEGER j3line(30)           !
      INTEGER j3dir(72)            ! connection direction used in only 1-D, + from 3-D to 1-D
      INTEGER j1Cvodn(72)          ! 3-D cell index at 1-D
      INTEGER j1Cvndx(72)          ! 3-D cell index at 1-D = i1Cvodn
      INTEGER j1jndx(72)           ! connection face index in 1-D face block
      INTEGER j1Rvodn(72)          ! interface 1-D cell index in cell block = i1Rvndx
      INTEGER j1Rvndx(72)          !
      INTEGER jndxsp(10,72)        !
      INTEGER j1nic(10,72)         !
      INTEGER s3sdbvolpack         !
      INTEGER s3sdbvolpacko(72)    !
      INTEGER s3c3pack(100,100)    !
      INTEGER s3vpackflag(100,100) !
      INTEGER s3hpackflag(100,100) !
      INTEGER s3overcorrection(72) !
      INTEGER j1max(10)
