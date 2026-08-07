!
      SUBROUTINE udfn_porous_region(x,ntype,ltype,xmin,xmax,property,arg,indx,nn) 
!
      USE Zparam        , ONLY: ndim
!
      IMPLICIT NONE      
! 
!
!     input
      INTEGER ltype,ntype,nn
      INTEGER indx(ltype)
      REAL(8) x(ndim)
      REAL(8) xmin(ndim,ltype),xmax(ndim,ltype)    
      REAL(8) property(nn,ntype)
!     output
      REAL(8) arg(nn)
!     local variables
      INTEGER ii
      REAL(8) xmin0(ndim),xmax0(ndim)
!      
      DO ii=1,ltype
           
          xmin0(:)=xmin(:,ii)  ! : means dimensions
          xmax0(:)=xmax(:,ii)
          IF(xmin0(1).lt.xmax0(1)) THEN  ! .lt. --> .and.
             IF(x(1).gt.xmin0(1).and.x(1).lt.xmax0(1)) THEN
                IF(xmin0(2).lt.xmax0(2)) THEN  ! .lt. --> .and.
                   IF(x(2).gt.xmin0(2).and.x(2).lt.xmax0(2)) THEN
!
                      IF(ndim.eq.3) THEN
                         IF(xmin0(3).lt.xmax0(3)) THEN  ! .lt. --> .and.
                            IF(x(3).gt.xmin0(3).and.x(3).lt.xmax0(3)) THEN                                 
                               !!!!!!!!!!!!!!
                               arg(:)=property(:,indx(ii))
                            ENDIF
!                                   
                         ELSEIF(xmin0(3).gt.xmax0(3)) THEN  ! .gt. --> .or.
                            IF(x(3).gt.xmin0(3).or.x(3).lt.xmax0(3)) THEN                                 
                               !!!!!!!!!!!!!!
                               arg(:)=property(:,indx(ii))
                            ENDIF
                         ENDIF
                      ELSE
                         !!!!!!!!!!!! 
                         arg(:)=property(:,indx(ii))
                      ENDIF   
!                              
                   ENDIF 
!                           
                ELSEIF(xmin0(2).gt.xmax0(2)) THEN ! .gt. --> .or.
                   IF(x(2).gt.xmin0(2).or.x(2).lt.xmax0(2)) THEN
!                              
                      IF(ndim.eq.3) THEN
                         IF(xmin0(3).lt.xmax0(3)) THEN  ! .lt. --> .and.
                            IF(x(3).gt.xmin0(3).and.x(3).lt.xmax0(3)) THEN                                 
                               !!!!!!!!!!!!!!
                               arg(:)=property(:,indx(ii))                                     
                            ENDIF
!                                   
                         ELSEIF(xmin0(3).gt.xmax0(3)) THEN  ! .gt. --> .or.
                            IF(x(3).gt.xmin0(3).or.x(3).lt.xmax0(3)) THEN                                 
                               !!!!!!!!!!!!!!
                               arg(:)=property(:,indx(ii))                                  
                            ENDIF
                         ENDIF
                      ELSE
                         !!!!!!!!!!!!!!  
                         arg(:)=property(:,indx(ii))                      
                      ENDIF                          
                        
                   ENDIF
                ENDIF
            ENDIF
!                     
         ELSEIF(xmin0(1).gt.xmax0(1)) THEN  ! .gt. --> .or.
            IF(x(1).gt.xmin0(1).or.x(1).lt.xmax0(1)) THEN
!
               IF(xmin0(2).lt.xmax0(2)) THEN  ! .lt. --> .and.
                  IF(x(2).gt.xmin0(2).and.x(2).lt.xmax0(2)) THEN
!
                     IF(ndim.eq.3) THEN
                        IF(xmin0(3).lt.xmax0(3)) THEN  ! .lt. --> .and.
                           IF(x(3).gt.xmin0(3).and.x(3).lt.xmax0(3)) THEN                                 
                              !!!!!!!!!!!!!!
                              arg(:)=property(:,indx(ii))                          
                           ENDIF
!                                   
                        ELSEIF(xmin0(3).gt.xmax0(3)) THEN  ! .gt. --> .or.
                           IF(x(3).gt.xmin0(3).or.x(3).lt.xmax0(3)) THEN                                 
                              !!!!!!!!!!!!!!
                              arg(:)=property(:,indx(ii))                                 
                           ENDIF
                        ENDIF
                     ELSE
                        !!!!!!!!!!!!!!   
                        arg(:)=property(:,indx(ii))                    
                     ENDIF   
!                              
                  ENDIF 
!                           
               ELSEIF(xmin0(2).gt.xmax0(2)) THEN ! .gt. --> .or.
                  IF(x(2).gt.xmin0(2).or.x(2).lt.xmax0(2)) THEN
!                              
                     IF(ndim.eq.3) THEN
                        IF(xmin0(3).lt.xmax0(3)) THEN  ! .lt. --> .and.
                           IF(x(3).gt.xmin0(3).and.x(3).lt.xmax0(3)) THEN                                 
                              !!!!!!!!!!!!!!
                              arg(:)=property(:,indx(ii))                                 
                           ENDIF
!                                   
                        ELSEIF(xmin0(3).gt.xmax0(3)) THEN  ! .gt. --> .or.
                           IF(x(3).gt.xmin0(3).or.x(3).lt.xmax0(3)) THEN                                 
                              !!!!!!!!!!!!!!
                              arg(:)=property(:,indx(ii))                                
                           ENDIF
                        ENDIF
                     ELSE
                        !!!!!!!!!!!!!!   
                        arg(:)=property(:,indx(ii))                          
                     ENDIF                          
!                       
                  ENDIF
               ENDIF
!                        
            ENDIF                     
         ENDIF
      ENDDO   
!      
      RETURN
      END SUBROUTINE udfn_porous_region      
