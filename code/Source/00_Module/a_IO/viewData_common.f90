      MODULE set_kind 
!     
      IMPLICIT NONE  
      SAVE 
!        
      INTEGER,PARAMETER :: dble = selected_real_kind(p=10,r=50)
      INTEGER,PARAMETER :: sngl = selected_real_kind(p=5)
      INTEGER,PARAMETER :: ddbl = selected_real_kind(p=10,r=50)
!        
      END MODULE set_kind 
!------------------------------------------------------------------------
      MODULE viewDataType
!
      USE set_kind
!      
      IMPLICIT NONE
      SAVE
!      
      INTEGER, PARAMETER :: MAXvPoint = 10, MAXvLine = 10, MAXvSurface = 10
      INTEGER, PARAMETER :: MAXvVar = 200
!      
      TYPE :: vectorType
         REAL(kind=dble) :: x,y,z
      END TYPE vectorType
!      
      TYPE :: lineType
         TYPE(vectorType) :: p1, p2
      END TYPE lineType  
!      
      TYPE :: viewPointType
         INTEGER:: nPoints
         TYPE(vectorType) :: coord(MAXvPoint)
         REAL(kind=dble) :: tplot, dtplot
         INTEGER:: nVariables
         CHARACTER(10):: var(MAXvVar)
      END TYPE viewPointType
!      
      TYPE :: viewLineType
         INTEGER :: nLines
         TYPE(lineType) :: line(MAXvLine)
         INTEGER :: nProbe
         REAL(kind=dble) :: tplot, dtplot
         INTEGER:: nVariables
         CHARACTER(10) :: var(MAXvVar)
      END TYPE viewLineType
!     
      TYPE :: viewSurfaceType
         INTEGER :: nSurfaces
         INTEGER :: surface(MAXvSurface)
         REAL(kind=dble) :: tplot, dtplot
         INTEGER :: nVariables
         CHARACTER(10) :: var(MAXvVar)
      END TYPE viewSurfaceType
!    
      TYPE :: viewFieldType0
         REAL(kind=dble) :: tplot, dtplot
         INTEGER :: nVariables
         CHARACTER(10) :: var(MAXvVar)
      END TYPE viewFieldType0
!      
      TYPE :: viewFieldType
         REAL(kind=dble) :: tplot, dtplot
         INTEGER :: nVectors, nScalars
         CHARACTER(20) :: vectorVar(MAXvVar), scalarVar(MAXvVar) 
      END TYPE viewFieldType
!      
      END MODULE viewDataType
!------------------------------------------------------------------------
      MODULE viewData_common 
!      
      USE set_kind
      USE viewDataType
!      
      TYPE(viewPointType) :: viewPoint
      TYPE(viewLineType) :: viewLine
      TYPE(viewSurfaceType) :: viewSurface
      TYPE(viewFieldType) :: viewField
!      
      CHARACTER(200) :: VFvariableNames
      CHARACTER(30) :: viwname
      INTEGER,PARAMETER :: viwUnit =100
      INTEGER :: nvector,nscalar
      REAL(8) :: crit_zero
      INTEGER,DIMENSION(:),ALLOCATABLE :: cupid_rv_jperm
      REAL(8),DIMENSION(:),ALLOCATABLE :: alphag_all,alphal_all,alphad_all
!      
      NAMELIST /postparam/ viewField
      REAL(kind=dble) :: time2plot
      INTEGER nframe
!      
      END MODULE viewData_common
