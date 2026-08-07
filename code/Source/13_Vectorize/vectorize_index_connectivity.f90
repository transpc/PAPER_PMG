!-------------------------------------------------------------------------------
      SUBROUTINE vectorize_index_connectivity
      USE Zvec_index       , ONLY: left_nf,right_fsw
      USE Zvec_index_solid , ONLY: left_solid_nf,right_solid_fsw,f_fluid_fsw,flux_fsw,qliq_fsw,qgas_fsw
      USE Zvec_index_solid , ONLY: dfilm_fsw,vfilm_fsw,dfilm_ctw,vfilm_ctw
      USE Zconst1          , ONLY: vv_prob
      USE Znum_cell        , ONLY: istart_nf,istartc_nf    
      USE Zvec_param       , ONLY: nf_fsw,nf_ctw
      USE Zcore            , ONLY: myrank
                                         
      IMPLICIT NONE
      INTEGER :: nf_number_c,istart_c,len_c,nf_number,istart,isize
      INTEGER :: i,i1,ii,kk,i_c,i1_c,ii_c,kk_c
!
!.....n_fluid_fsw
!
      IF(nf_fsw.le.0)GOTO 1 
      nf_number_c=1
      istart_c=istartc_nf(1,nf_number_c)
      len_c   =istartc_nf(2,nf_number_c)
!
      nf_number=5
      istart=istart_nf(1,nf_number)
      isize =istart_nf(2,nf_number)
!
      IF(len_c.ne.nf_fsw)THEN
         WRITE(*,"(11x,a,3i5)")'Error in len_c,nf_fsw,myrank:',len_c,nf_fsw,myrank
         STOP
      ENDIF  
      IF(isize.ne.nf_fsw)THEN
         WRITE(*,"(11x,a,3i5)")'Error in isize,nf_fsw,myrank:',isize,nf_fsw,myrank
         STOP
      ENDIF   
      ALLOCATE(f_fluid_fsw(nf_fsw))
      ALLOCATE(flux_fsw(nf_fsw))
      ALLOCATE(qliq_fsw(nf_fsw))
      ALLOCATE(qgas_fsw(nf_fsw))
      f_fluid_fsw(:)=0
      flux_fsw(:)=0.0d0
      qliq_fsw(:)=0.d0
      qgas_fsw(:)=0.d0         
      !ALLOCATE(tliq_fsw(nf_fsw))
      !ALLOCATE(tsol_fsw(nf_fsw))      
      !tliq_fsw(:)=0.0d0
      !tsol_fsw(:)=0.0d0
! 
      IF(vv_prob.eq.'kuhn_111'.or.vv_prob.eq.'Nuscale-03Pool')THEN        
         ALLOCATE(dfilm_fsw(nf_fsw))
         ALLOCATE(vfilm_fsw(nf_fsw))
         dfilm_fsw(:)=0.0d0
         vfilm_fsw(:)=0.0d0   
      ENDIF      
!            
      DO i=1,isize
         i1=istart+i
         ii=left_nf(i1)
         kk=right_fsw(i)
!         
         DO i_c=1,len_c
            i1_c=istart_c+i_c
            ii_c=left_solid_nf(i1_c)
            kk_c=right_solid_fsw(i_c)     !Right: Fluid Side          
            IF(kk_c.eq.ii)THEN
            	 IF(f_fluid_fsw(i_c).eq.0)THEN
                  f_fluid_fsw(i_c)=i
                  EXIT
               ENDIF
            ENDIF
         ENDDO      
!            
      ENDDO         
!      
1     CONTINUE
!
      IF((vv_prob.eq.'kuhn_111'.or.vv_prob.eq.'Nuscale-03Pool').and.nf_ctw.gt.0)THEN
         ALLOCATE(dfilm_ctw(nf_ctw))
         ALLOCATE(vfilm_ctw(nf_ctw))      
         dfilm_ctw(:)=0.0d0
         vfilm_ctw(:)=0.0d0
      ENDIF
!      
      RETURN
!      
      END SUBROUTINE vectorize_index_connectivity      
