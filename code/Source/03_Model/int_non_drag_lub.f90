!
      SUBROUTINE int_non_drag_lub
!
!     This routine calculates wall lubrication force
!
      USE VOL_DATA                    
      USE Zzone           , ONLY: ncell_fluid,ncell_fluid_all
      USE Zcore           , ONLY: np      
      USE Zparam          , ONLY: ndim
      USE Zbc_index       , ONLY: num_wall_group
      USE Zndforce        , ONLY: cwlf,d_bfc,F_wl,face_wall_group,dis_closewall,cell_closewall, &
                                  cell_closewall_indx
      USE Zndforce        , ONLY: c_bface_indx,s_wall_lub
      USE Znormal         , ONLY: xn_wallcell
      USE Zvector         , ONLY: vg_n,vl_n
      USE Zuserdefined    , ONLY: udfl_psbt_cfx_model
      
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i,k,ii,jj,m,j0
      REAL(8) :: d_bubble,Awlf,vrela(ndim),vrela2
      REAL(8) :: Cwlf1,Cwlf2,Cwlf3,d_wall
      REAL(8) :: factor
!.....Local arrays
      REAL(8) :: Awlf_v(ncell_fluid)
      INTEGER, ALLOCATABLE :: face_wall_group_all(:)
!      
!.....Set the model coefficients
!     
      IF(udfl_psbt_cfx_model)THEN
         Cwlf1=-0.025d0
         Cwlf2=0.075d0
         Cwlf3=0.0d0
      ELSE
          IF(s_wall_lub.eq.'antal&cfx')THEN                      !cfx
             Cwlf1=-0.01d0
             Cwlf2=+0.05d0 
             Cwlf3=0.0d0                                        
          ELSEIF(s_wall_lub.eq.'antal&antal')THEN               !antal original
             Cwlf1=-0.104d0
             Cwlf2=+0.147d0   
             Cwlf3=-0.06                                  
          ENDIF 
      ENDIF       
! 
      IF(s_wall_lub.eq.'antal&antal') THEN
         factor=0.5d0
      ELSE
         factor=1.d0
      ENDIF
!
      F_wl(:,:)=0.0d0
!
!.....If num_wall_group=0, use the d_bfc which is the most close wall distance
!           
      IF(num_wall_group.eq.0)THEN
!          
!........Antal model
!
         IF(s_wall_lub.eq.'antal&cfx'.or.s_wall_lub.eq.'antal&antal')THEN 
            DO i=1,ncell_fluid
               j0=c_bface_indx(i)
               d_bubble=cell%d1(i)
               d_bubble=d_bubble*factor
               vrela(:)=(vg_n(i,:)-vl_n(i,:))-DOT_PRODUCT(xn_wallcell(j0,:),(vg_n(i,:)-vl_n(i,:)))*xn_wallcell(j0,:)
               vrela2=DOT_PRODUCT(vrela(:),vrela(:))
               IF(d_bfc(i).eq.0)THEN
                  Awlf=0.0d0
               ELSE
                  Awlf=cell%alphag(i)*cell%rhol(i)*Cwlf(i)*vrela2/d_bubble * &
                       MAX(0.d0,Cwlf1+Cwlf2*d_bubble/d_bfc(i)+Cwlf3*SQRT(vrela2))            
               ENDIF
               F_wl(i,1)=Awlf*xn_wallcell(j0,1)                                   ! Consider only one wall face 
               F_wl(i,2)=Awlf*xn_wallcell(j0,2)  
               IF(ndim.eq.3) F_wl(i,3)=Awlf*xn_wallcell(j0,3) 
            ENDDO   
!             
!........Tomiyama model    
!
         ELSEIF(s_wall_lub.eq.'tomiyama')THEN 
            CALL udfn_wall_lub_v(d_bfc,Awlf_v)  
            DO i=1,ncell_fluid
               j0=c_bface_indx(i)
               d_bubble=cell%d1(i)
               d_bubble=d_bubble*factor
               vrela(:)=(vg_n(i,:)-vl_n(i,:))-DOT_PRODUCT(xn_wallcell(j0,:),(vg_n(i,:)-vl_n(i,:)))*xn_wallcell(j0,:)
               vrela2=DOT_PRODUCT(vrela(:),vrela(:))
               Awlf=Awlf_v(i) 
               F_wl(i,1)=Awlf*xn_wallcell(j0,1)                                   ! Consider only one wall face 
               F_wl(i,2)=Awlf*xn_wallcell(j0,2)  
               IF(ndim.eq.3) F_wl(i,3)=Awlf*xn_wallcell(j0,3) 
            ENDDO   
         ENDIF
!          
      ELSE
!
!.....If num_wall_group>0, use dis_closewall(i,m) which is the most close wall distances to the each wall group
!        
         ALLOCATE(face_wall_group_all(ncell_fluid_all)) 
         DO m=1,num_wall_group
!
            face_wall_group_all(:)=0 
            IF (np.gt.1) THEN
               DO k=1,num_wall_group
                  CALL allgatherv_i(face_wall_group(1,k),face_wall_group_all,ncell_fluid,ncell_fluid_all,0)
               ENDDO
            ELSE
               face_wall_group_all(:)=face_wall_group(:,k)
            ENDIF      
!          
!...........Antal model
!
            IF(s_wall_lub.eq.'antal&cfx'.or.s_wall_lub.eq.'antal&antal')THEN 
               DO i=1,ncell_fluid
                  d_wall=dis_closewall(i,m)
                  ii=cell_closewall(i,m)                                          ! ii is grobal cell number
                  jj=face_wall_group_all(ii)
                  d_bubble=cell%d1(i)
                  IF(jj.eq.0.or.ii.eq.0.or.d_wall.eq.0)THEN
                  ELSE
                     j0=cell_closewall_indx(i,m)
                     d_bubble=d_bubble*factor
                     vrela(:)=(vg_n(i,:)-vl_n(i,:))-DOT_PRODUCT(xn_wallcell(j0,:),(vg_n(i,:)-vl_n(i,:)))*xn_wallcell(j0,:)
                     vrela2=DOT_PRODUCT(vrela(:),vrela(:))
                     Awlf=cell%alphag(i)*cell%rhol(i)*Cwlf(i)*vrela2/d_bubble * &
                          MAX(0.d0,Cwlf1+Cwlf2*d_bubble/d_wall+Cwlf3*SQRT(vrela2))
                     F_wl(i,1)=F_wl(i,1)+Awlf*xn_wallcell(j0,1)  
                     F_wl(i,2)=F_wl(i,2)+Awlf*xn_wallcell(j0,2)  
                     IF(ndim.eq.3) F_wl(i,3)=F_wl(i,3)+Awlf*xn_wallcell(j0,3)
                  ENDIF
               ENDDO    
!             
!...........Tomiyama model    
!
            ELSEIF(s_wall_lub.eq.'tomiyama')THEN
               DO i=1,ncell_fluid
                  d_wall=dis_closewall(i,m)
                  ii=cell_closewall(i,m)                                          ! ii is grobal cell number
                  jj=face_wall_group_all(ii)
                  d_bubble=cell%d1(i)
                  IF(jj.eq.0.or.ii.eq.0.or.d_wall.eq.0)THEN
                  ELSE
                     j0=cell_closewall_indx(i,m)
                     d_bubble=d_bubble*factor
                     vrela(:)=(vg_n(i,:)-vl_n(i,:))-DOT_PRODUCT(xn_wallcell(j0,:),(vg_n(i,:)-vl_n(i,:)))*xn_wallcell(j0,:)
                     vrela2=DOT_PRODUCT(vrela(:),vrela(:))
                     CALL udfn_wall_lub_i(i,d_wall,Awlf)    
                     F_wl(i,1)=F_wl(i,1)+Awlf*xn_wallcell(j0,1)  
                     F_wl(i,2)=F_wl(i,2)+Awlf*xn_wallcell(j0,2)  
                     IF(ndim.eq.3) F_wl(i,3)=F_wl(i,3)+Awlf*xn_wallcell(j0,3)
                  ENDIF
               ENDDO    
            ENDIF
!
         ENDDO
         DEALLOCATE(face_wall_group_all)         
      ENDIF           
!        
      END SUBROUTINE int_non_drag_lub
