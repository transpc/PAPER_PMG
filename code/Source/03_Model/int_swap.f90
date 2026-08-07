!
      SUBROUTINE int_swap(iflag)
!
      USE Vol_DATA       , ONLY: cell
      USE Zzone          , ONLY: ncell_fluid,icore
      USE Zqvol          , ONLY: H_ig,H_il,H_gf,hil_o,hig_o
      USE Zrv_model      , ONLY: rv_model,                        &
                                 lfric_swap,lhtc_swap,lfricw_swap
      USE Zmodel         , ONLY: vfgl_cfd,vfgd_cfd,h_il_cfd,h_ig_cfd,h_gf_cfd,cd_min_ag99,cd_min_user,cb_bubble
      USE Zndforce       , ONLY: relax_hik,relax_cd,vfgl_o
      USE Zuserdefined   , ONLY: udfl_model_overwrite
      
!      
      IMPLICIT NONE
!
!.....Input
      INTEGER iflag
!.....Local variables
      INTEGER i
!
!.....Initialize variables or store interfacial drag and heat transfer coefficients of open media.
!      
      IF(iflag.eq.0)THEN
         Cb_bubble(:)=0.44d0
         cell%entr(:)=0.0d0
         cell%dentr(:)=0.0d0
         cell%yeta(:)=0.0d0  
         RETURN
      ELSEIF(iflag.eq.1.and.rv_model.gt.0)THEN
         vfgl_cfd(:)=cell%vfgl(:)
         vfgd_cfd(:)=cell%vfgd(:)
         h_il_cfd(:)=h_il(:)
         h_ig_cfd(:)=h_ig(:)      
         h_gf_cfd(:)=h_gf(:)   
         RETURN   
      ELSEIF(iflag.eq.2.and.rv_model.gt.0)THEN
         h_il_cfd(:)=h_il(:)
         h_ig_cfd(:)=h_ig(:)      
         h_gf_cfd(:)=h_gf(:)   
         RETURN              
      ENDIF
!
!.....Interfacial heat transfer coefficient
!
      IF(lhtc_swap.and.rv_model.gt.0)THEN
         DO i=1,ncell_fluid 
            IF(icore(i).ne.1)THEN
                h_il(i)=h_il_cfd(i)
                h_ig(i)=h_ig_cfd(i)  
                h_gf(i)=h_gf_cfd(i)
            ENDIF    
         ENDDO
      ENDIF
!.....relaxation   
      IF(relax_hik.gt.1.0d-10)THEN      
         H_il(:)=(1.0d0-relax_hik)*H_il(:)+relax_hik*hil_o(:)
         H_ig(:)=(1.0d0-relax_hik)*H_ig(:)+relax_hik*hig_o(:)
         !IF(free_model.and.mhtc.eq.0)THEN !next
         !   !See int_htc_simple_model.f90.
         !ELSE
         !   !See int_htc_full_topology.f90,int_htc_simple_topology.f90
         !   H_gf(:)=H_ig(:)
         !ENDIF   
      ENDIF         
!      
      IF(iflag.ne.11)RETURN      
! 
!.....Interfacial drag coefficient
!
     IF(lfric_swap.and.rv_model.gt.0)THEN  
        DO i=1,ncell_fluid 
            IF(icore(i).ne.1)THEN
                cell%vfgl(i)=vfgl_cfd(i)
                cell%vfgd(i)=vfgd_cfd(i)  
            ENDIF    
         ENDDO
     ENDIF
!     
!....relaxation      
!     
     IF(relax_cd.eq.0.0d0)THEN
        DO i=1,ncell_fluid
           cell%vfgl(i)=MAX(cell%vfgl(i),Cd_min_user)
           IF(cell%alphag(i).gt.0.99d0)cell%vfgl(i)=MAX(cell%vfgl(i),Cd_min_ag99)
           IF(cell%alphag(i).lt.0.01d0)cell%vfgl(i)=MAX(cell%vfgl(i),Cd_min_ag99)
        ENDDO        
     ELSE
        DO i=1,ncell_fluid
           cell%vfgl(i)=vfgl_o(i)*relax_cd+cell%vfgl(i)*(1.d0-relax_cd)
           cell%vfgl(i)=MAX(cell%vfgl(i),Cd_min_user)
           IF(cell%alphag(i).gt.0.99d0)cell%vfgl(i)=MAX(cell%vfgl(i),Cd_min_ag99)
           IF(cell%alphag(i).lt.0.01d0)cell%vfgl(i)=MAX(cell%vfgl(i),Cd_min_ag99)
        ENDDO          
     ENDIF
! 
!.....Wall friction coefficient
!
!      IF(rv_fric_w.gt.0.and.udfl_icore)THEN
      IF(lfricw_swap)THEN
         DO i=1,ncell_fluid 
            IF(icore(i).ne.1)THEN
                cell%vfwl(i)=0.0d0
                cell%vfwg(i)=0.0d0   
            ENDIF    
         ENDDO
      ENDIF   
!
!.....Overwrite the coefficients of the models
!
      IF(udfl_model_overwrite)CALL udfn_model_overwrite
!            
      END SUBROUTINE int_swap

