!
      SUBROUTINE cupid2mars 
!
!     This routines transfer informations from CUPID to MARS,
!      
      USE Vol_DATA
      USE Zzone      , only:ncell_fluid
      USE Zcore      , only:np,myrank
      USE Ztimecon   , ONLY:time
      USE Zmars      , only:marsindex,n_marsbc,pcup,egcup,elcup,&
                            alphagcup,rholcup,rhogcup,cboroncup,qualacup,&
                            pcup_tmp,egcup_tmp,elcup_tmp,alphagcup_tmp,&
                            rhogcup_tmp,rholcup_tmp,cboroncup_tmp,qualacup_tmp,&
                            tlcup,tgcup,tlcup_tmp,tgcup_tmp
      USE Zmpi       , only:jperm
      USE c3com_cupid, only:i3cupid

!
      IMPLICIT none
!      
!DEC$IF defined (MCC)      
!DEC$ELSEIF defined (MCC_DLL)      
      !dec$ attributes dllexport :: cupid2mars      
!DEC$ELSEIF defined (SPACE)          
      !dec$ attributes dllexport :: cupid2mars
!DEC$ENDIF
!  
      INCLUDE 'c3com.h' 
!
      INTEGER i3ii,i3jj,i,j
!
      DO i=1,ncell_fluid
         IF(marsindex(i).ge.1 .and. marsindex(i).le.n_marsbc)then
!             
            IF(i3cupid(marsindex(i)).ne.jperm(i))then !marsindex(i)=k
               WRITE(*,*)'global cell number is not coherent!', i3cupid(marsindex(i)),jperm(i)
               STOP
            ENDIF   
!            
            pcup(marsindex(i))=cell%p(i)
            egcup(marsindex(i))=cell%eg(i)
            elcup(marsindex(i))=cell%el(i)
            alphagcup(marsindex(i))=cell%alphag(i)
            rholcup(marsindex(i))=cell%rhol(i)
            rhogcup(marsindex(i))=cell%rhog(i)
            cboroncup(marsindex(i))=cell%cboron(i)
            qualacup(marsindex(i))=cell%quala(i)
            tgcup(marsindex(i))=cell%tg(i)
            tlcup(marsindex(i))=cell%tl(i)
!            
         ENDIF     
      ENDDO
!      
      CALL allreduce_r(pcup     ,pcup_tmp     ,n_marsbc)
      CALL allreduce_r(egcup    ,egcup_tmp    ,n_marsbc)
      CALL allreduce_r(elcup    ,elcup_tmp    ,n_marsbc)
      CALL allreduce_r(alphagcup,alphagcup_tmp,n_marsbc)
      CALL allreduce_r(rholcup  ,rholcup_tmp  ,n_marsbc)
      CALL allreduce_r(rhogcup  ,rhogcup_tmp  ,n_marsbc)
      CALL allreduce_r(cboroncup,cboroncup_tmp,n_marsbc)
      CALL allreduce_r(qualacup ,qualacup_tmp ,n_marsbc)
      CALL allreduce_r(tgcup    ,tgcup_tmp    ,n_marsbc)
      CALL allreduce_r(tlcup    ,tlcup_tmp    ,n_marsbc)
!                                                                       
!...Data transfer from CUPID to MARS                 
!                                                                       
      i3ii=1 
      DO i3jj=1,i3nic(2) 
         j=i3cupid(i3jj) 
!                                                                    
!......pressure [Pa]                                      
!                                                                    
         c3pa(i3ii,i3jj)=pcup_tmp(i3jj) 
!                                                                    
!......internal energy [J/kg]                                  
!
         c3ug(i3ii,i3jj)=egcup_tmp(i3jj) 
         c3uf(i3ii,i3jj)=elcup_tmp(i3jj) 
!
!......void fraction [-]                                                  
!                                                                    
         c3al(i3ii,i3jj)=alphagcup_tmp(i3jj) 
! 
!......density [kg/m3]                
!                                                                    
         c3rhof(i3ii,i3jj)=rholcup_tmp(i3jj) 
         c3rhog(i3ii,i3jj)=rhogcup_tmp(i3jj)
!
!......boron concentration[1] 
!
         c3brn(i3ii,i3jj)=cboroncup_tmp(i3jj)
!                                                                    
!......noncondensible gas mass / volume [kg/m3]                
!                                                                    
         c3arxq(i3ii,i3jj)=qualacup_tmp(i3jj) 
!
!......sdbvol vapor generation [kg/s.m3]
!
         c3vpgno(i3ii,i3jj)=0.0d0   !MCC-jjj-NEXT
!         
      ENDDO 
!
!DEC$IF defined (SPACE)
      CALL c3com_copy_C2S
!DEC$ENDIF
      RETURN 
      END SUBROUTINE cupid2mars                            
