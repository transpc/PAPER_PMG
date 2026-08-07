      MODULE STM_TBL_cupid
      IMPLICIT NONE
      SAVE
!
!     NC gases are: Helium, Hydrogen, Nitrogen, Krypton, Xenon, Air,
!     Argon, SF6.
      REAL(8) :: wmole(8),visa(8),thcax(8),thcbx(8),dcvax(8),cvaox(8),uaox(8),advn(8)
      DATA wmole/4.002598d0,2.01593d0,28.01403d0,83.800d0,131.300d0,28.963d0,39.948d0,146.05d0/
      DATA visa/1.473d-6,6.675d-7,1.381d-6,2.386d-6,3.455d-6,1.492d-6,1.935d-6,2.306654d-6/
      DATA thcax/2.639d-3,1.097d-3,5.314d-4,8.247d-5,4.351d-5,1.945d-4,2.986d-4,2.374568d-2/
      DATA thcbx/0.7085d0,0.8785d0,0.6898d0,0.8363d0,0.8616d0,0.8586d0,0.7224d0,0.0d0/
      DATA dcvax/0.003455924d0,0.522573d0,0.1184518d0,0.0035d0,0.0035d0,0.10329037d0,0.003517d0,1.0d-6/
      DATA cvaox/3115.839d0,10310.75d0,741.9764d0,148.824d0,94.9084d0,715.0d0,312.192d0,793.399d0/
      DATA uaox/13256.44d0,182783.4d0,145725.884d0,122666.5d0,122666.5d0,158990.52d0,122666.5d0,0.0d0/
      DATA advn/2.67d0,6.12d0,18.5d0,24.5d0,32.7d0,19.7d0,16.2d0,71.3d0/ 
!
      REAL(8) :: rmolg=8314.3d0
!
      REAL(8) :: pxxx=1378.951459d0
      REAL(8) :: pxxy=3102640.782d0
      REAL(8) :: pxx1=1.450377377d-3
      REAL(8) :: pxx2=1.450377377d-4
!
      REAL(8) :: crt=647.3d0
      REAL(8) :: crp=22120000.d0
!
      REAL(8) :: b(6),c(9),cc(3),k(9),g(5)
      DATA b/6669.352222d0,-4658.899d0,1376.536722d0,-201.9126167d0, &
             14.82832111d0,-.4337434056d0/
      DATA c/274.9043833d0,13.66254889d0,1.176781611d0,-.189693d0,       &
             8.74535666d-2,-1.7405325d-2,2.147682333d-3,-1.383432444d-4, &
             3.800086611d-6/
      DATA cc/0.84488898d0,2.9056480d0,219.74589d0/
      DATA k/-7.691234564d0,-26.08023696d0,-168.1706546d0,6.423285504d1, &
             -1.189646225d2,4.167117320d0,2.097506760d1,1.d9,6.0d0/
!.....g(i)=i*k(i) i=1,5
      DATA g/-7.691234564d0,-52.16047392d0,-504.5119638d0,25.693142016d1, &
             -5.948231125d2/
!      
!.....steam table block
!
      REAL(8),DIMENSION(:),ALLOCATABLE :: st_tbl
      REAL(8),DIMENSION(:),ALLOCATABLE :: a31,a41
      REAL(8),DIMENSION(:,:,:),ALLOCATABLE :: a3,a4
!
!.....Old sth2xc common block set in read_tpfh2o
!      
      INTEGER :: nt,np,ns,ns2
!
!.....Old stcblk common block ndxstd set in read_tpfh2o
!     nfluid set in read_flow either default=1 or in somaFlow
!
      INTEGER :: ndxstd
!
!  nfluid  fluid number:
!           1 = light water
!           2 = heavy water
!           3 = carbon dioxide
!           4 = helium
!           5 = hydrogen
!           6 = oxygen
!           7 = nitrogen
!           8 = sodium
!           9 = lead bismuth eutetic
!          10 = reserved
!          11 = reserved
!          12 = new light water (based on internal energy and pressure)
!          13 = reserved
!          14 = reserved
!          15 = new liqht water (supercritical etc. - restructured)
!          16 = r12
!          17 = r134a
      INTEGER :: nfluid
!
!.....Old stcom
!
      REAL (8) :: pcrit
!
!.....Old NONCONG
!
      REAL(8) :: tao,cvao,uao,dcva,ra
!      
      END MODULE STM_TBL_cupid
