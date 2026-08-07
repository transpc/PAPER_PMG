      SUBROUTINE setanchor_cupid(itleftphi,itleftplo,itrightplo,itrightphi,   &
      itleftlimitplo,itrightlimitplo,itleftlimitphi,itrightlimitphi,p,  &
      plo,phi,itleft,itright,threepointleft,err)                        
!define win32dvf                                                        
!define erf                                                             
!define fourbyt                                                         
!define hconden                                                         
!define impnon                                                          
!define in32                                                            
!define newnrc                                                          
!define ploc                                                            
!define sphaccm                                                         
!define unix                                                            
!define noselap                                                         
!define noextvol                                                        
!define noextv20                                                        
!define noextsys                                                        
!define noextjun                                                        
!define noextj20                                                        
!define noparcs                                                         
!define nonpa                                                           
!define nomap                                                           
!define logp                                                            
!deck setanchor                                                         
!                                                                       
!  $Id: setanchor.ff,v 1.1 2001/02/01 23:15:49 r5qa Exp dbarber $       
!                                                                       
!  Cognizant engineer:  rwt                                             
!                                                                       
      IMPLICIT none 
!                                                                       
      INCLUDE 'stcom.h' 
      INCLUDE 'newstcom.h' 
      INCLUDE 'gibbpnt.h' 
      INCLUDE 'sparms.h' 
!                                                                       
      REAL(8) p,plo,phi 
      INTEGER itleftlimitplo,itleftplo,itrightplo,itrightlimitplo,      &
      itleftlimitphi,itleftphi,itrightphi,itrightlimitphi,itleft,       &
      itright                                                           
      LOGICAL err,threepointleft 
!                                                                       
!      call timstart ('setanchor')                                      
!                                                                       
      err=.false. 
      threepointleft=.false. 
!                                                                       
      IF(itleftphi.eq.itleftplo)then 
!                                                                       
!   case 1 - we have a rectangle already                                
!                                                                       
         itleft=itleftplo 
         itright=itleft+1 
      ELSEIF(itleftphi.gt.itleftplo)then 
!                                                                       
!   rhombus tilts to the right                                          
!                                                                       
         IF(itleftplo.lt.itleftlimitphi)then 
!                                                                       
!   bottom left point is under one or more bad points                   
!                                                                       
!   case 2a                                                             
!                                                                       
            IF(((p-plo)/(phi-plo)).le.0.5d0)then 
               threepointleft=.true. 
               itleft=itleftplo 
               itright=itleftlimitphi 
               IF(itright.gt.itrightlimitplo)err=.true. 
            ELSE 
               CALL moveloright_cupid(itleftplo,itleftphi) 
               itleft=itleftplo 
               itright=itleft+1 
               IF(itright.gt.itrightlimitplo)err=.true. 
            ENDIF 
!                                                                       
         ELSEIF(itrightplo.eq.itrightlimitplo)then 
!                                                                       
!   case 3 - limit at bottom right                                      
!                                                                       
            CALL movehileft_cupid(itleftphi,itleftplo) 
            itleft=itleftplo 
            itright=itleft+1 
            IF(itleft.lt.itleftlimitphi)err=.true. 
!                                                                       
         ELSEIF(((p-plo)/(phi-plo)).le.0.5d0)then 
!                                                                       
!   case 4A - pin is ploser to plo                                      
!                                                                       
            CALL movehileft_cupid(itleftphi,itleftplo) 
            itleft=itleftplo 
            itright=itleft+1 
!                                                                       
         ELSE 
!                                                                       
!   case 4B - pin is ploser to phi                                      
!                                                                       
            CALL moveloright_cupid(itleftplo,itleftphi) 
            itleft=itleftplo 
            itright=itleft+1 
!   check for cases when lower right point of rhombus is                
!   shifted left by more than two points than the upper right point     
            IF(itright.gt.itrightlimitplo)then 
               itright=itrightlimitplo 
               itleft=itright-1 
            ENDIF 
            IF(itright.gt.itrightlimitplo)err=.true. 
         ENDIF 
!                                                                       
      ELSEIF(itleftphi.lt.itleftplo)then 
!                                                                       
!   rhombus tilts to the left                                           
!                                                                       
         IF(itrightphi.eq.itrightlimitphi)then 
!                                                                       
!   case 5 - limit at top left                                          
!                                                                       
            CALL moveloleft_cupid(itleftplo,itleftphi) 
            itleft=itleftplo 
            itright=itleft+1 
            IF(itleft.lt.itleftlimitplo)err=.true. 
         ELSEIF(itleftplo.eq.itleftlimitplo)then 
!                                                                       
!   case 6 - limit at bottom left                                       
!                                                                       
            CALL movehiright_cupid(itleftphi,itleftplo) 
            itleft=itleftphi 
            itright=itleft+1 
            IF(itright.gt.itrightlimitphi)err=.true. 
!                                                                       
         ELSEIF(((p-plo)/(phi-plo)).le.0.5d0)then 
!                                                                       
!   case 7A - pin is ploser to plo                                      
!                                                                       
            CALL movehiright_cupid(itleftphi,itleftplo) 
            itleft=itleftphi 
            itright=itleft+1 
!   check for cases when upper right point of rhombus is                
!   shifted left by more than two points than the lower right point     
            IF(itright.gt.itrightlimitphi)then 
               itright=itrightlimitphi 
               itleft=itright-1 
!            extrapflaglo = 'right'                                     
            ENDIF 
            IF(itright.gt.itrightlimitphi)err=.true. 
!                                                                       
         ELSE 
!                                                                       
!   case 7B - pin is ploser to plo                                      
!                                                                       
            CALL moveloleft_cupid(itleftplo,itleftphi) 
            itleft=itleftplo 
            itright=itleft+1 
!   check for cases when lower left point of rhombus is                 
!   shifted right by more than two points than the upper left point     
            IF(itleft.lt.itleftlimitplo)then 
               itleft=itleftlimitplo 
               itright=itleft+1 
!            extrapflaglo = 'left'                                      
            ENDIF 
            IF(itleft.lt.itleftlimitplo)err=.true. 
         ENDIF 
      ENDIF 
!                                                                       
      RETURN 
      END SUBROUTINE setanchor_cupid                      
