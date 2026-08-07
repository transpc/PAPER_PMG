      SUBROUTINE grad_limiter(s,dsdx,ilim)
!
!     This routine calculates the slope limiters for second-order
!     treatment of convection terms
!
      USE Zmpi         , ONLY: ncell_fp
      USE Zzone        , ONLY: ncell_fluid
      USE Zparam       , ONLY: ndim
      USE Z2nd_order   , ONLY: faclim,vkt
      USE Zcoord3      , ONLY: vol
      USE Znum_cell    , ONLY: istart_nf
      USE Zvec_index   , ONLY: left_nf,right_non
      USE Zvec_geo     , ONLY: dxfc_nf,    &
                               dxfc_non_k
!
      IMPLICIT NONE
!     input
      INTEGER ilim
      REAL(8) s(ncell_fp)
!     output
      REAL(8) dsdx(ncell_fp,ndim)
!     local variables
      INTEGER i,ii,kk,i1
      INTEGER nf_number,len,istart
      REAL(8) dels,dpp,dmm,psi,fs
      REAL(8) delh,eps2,tmp
!     local arrays
      REAL(8) spsi(ncell_fluid),smin(ncell_fluid),smax(ncell_fluid)
!
!.... Initialize limiter
!
      DO i=1,ncell_fluid
         spsi(i)=1.d0
         smin(i)=s(i)
         smax(i)=s(i)
      ENDDO

!
!.....Obtain min and max value among neighboring cells
!
      nf_number=0
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      DO i=1,len  
         i1=istart+i
         ii=left_nf(i1)
         kk=right_non(i)
         smax(ii)=DMAX1(smax(ii),s(kk))
         smin(ii)=DMIN1(smin(ii),s(kk))
         IF(kk.le.ncell_fluid) THEN
         smax(kk)=DMAX1(smax(kk),s(ii))
         smin(kk)=DMIN1(smin(kk),s(ii))
         ENDIF
      ENDDO
!
!.....Barth's limiter : AIAA Paper 89-0366
!
      IF(ilim.eq.1)THEN
!
         IF(ndim.eq.2) THEN
            nf_number=0
            len   =istart_nf(2,nf_number)
            DO i=1,len  
               i1=istart+i
               ii=left_nf(i1)
               kk=right_non(i)
!
               dels= dsdx(ii,1)*dxfc_nf(i1,1) &
                    +dsdx(ii,2)*dxfc_nf(i1,2)
               dpp=smax(ii)-s(ii)
               dmm=smin(ii)-s(ii)
               CALL barth_limiter(dels,dpp,dmm,psi)
               spsi(ii)=DMIN1(psi,spsi(ii))
!
         IF(kk.le.ncell_fluid) THEN
               dels= dsdx(kk,1)*dxfc_non_k(i,1) &
                    +dsdx(kk,2)*dxfc_non_k(i,2)
               dpp=smax(kk)-s(kk)
               dmm=smin(kk)-s(kk)
               CALL barth_limiter(dels,dpp,dmm,psi)
               spsi(kk)=DMIN1(psi,spsi(kk))
         ENDIF
            ENDDO
!
!........The rest
            DO nf_number=2,8
               istart=istart_nf(1,nf_number)
               len   =istart_nf(2,nf_number)
               DO i=1,len  
                  i1=istart+i
                  ii=left_nf(i1)
                  dels= dsdx(ii,1)*dxfc_nf(i1,1) &
                       +dsdx(ii,2)*dxfc_nf(i1,2)
                  dpp=smax(ii)-s(ii)
                  dmm=smin(ii)-s(ii)
                  CALL barth_limiter(dels,dpp,dmm,psi)
                  spsi(ii)=DMIN1(psi,spsi(ii))
               ENDDO
            ENDDO
         ELSE
            nf_number=0
            len   =istart_nf(2,nf_number)
            DO i=1,len  
               i1=istart+i
               ii=left_nf(i1)
               kk=right_non(i)
!
               dels= dsdx(ii,1)*dxfc_nf(i1,1) &
                    +dsdx(ii,2)*dxfc_nf(i1,2) &
                    +dsdx(ii,3)*dxfc_nf(i1,3)
               dpp=smax(ii)-s(ii)
               dmm=smin(ii)-s(ii)
               CALL barth_limiter(dels,dpp,dmm,psi)
               spsi(ii)=DMIN1(psi,spsi(ii))
!
         IF(kk.le.ncell_fluid) THEN
               dels= dsdx(kk,1)*dxfc_non_k(i,1) &
                    +dsdx(kk,2)*dxfc_non_k(i,2) &
                    +dsdx(kk,3)*dxfc_non_k(i,3)
               dpp=smax(kk)-s(kk)
               dmm=smin(kk)-s(kk)
               CALL barth_limiter(dels,dpp,dmm,psi)
               spsi(kk)=DMIN1(psi,spsi(kk))
         ENDIF
            ENDDO
!
!........The rest
            DO nf_number=2,8
               istart=istart_nf(1,nf_number)
               len   =istart_nf(2,nf_number)
               DO i=1,len  
                  i1=istart+i
                  ii=left_nf(i1)
                  dels= dsdx(ii,1)*dxfc_nf(i1,1) &
                       +dsdx(ii,2)*dxfc_nf(i1,2) &
                       +dsdx(ii,3)*dxfc_nf(i1,3)
                  dpp=smax(ii)-s(ii)
                  dmm=smin(ii)-s(ii)
                  CALL barth_limiter(dels,dpp,dmm,psi)
                  spsi(ii)=DMIN1(psi,spsi(ii))
               ENDDO
            ENDDO
         ENDIF
!
!.....Venkatakrishnan's limiter : JCP Paper 93-0880
!
      ELSEIF(ilim.eq.2)THEN
!
         IF(ndim.eq.2) THEN
            nf_number=0
            len   =istart_nf(2,nf_number)
            DO i=1,len  
               i1=istart+i
               ii=left_nf(i1)
               kk=right_non(i)
!
               delh=vol(ii)**(1.d0/2.d0)
               IF(vkt.gt.0.0d0)THEN
                  tmp=vkt*delh
                  eps2=tmp*tmp*tmp
               ELSE
                 eps2=0.0d0
               ENDIF
!
               dels= dsdx(ii,1)*dxfc_nf(i1,1) &
                    +dsdx(ii,2)*dxfc_nf(i1,2)
               dpp=smax(ii)-s(ii)
               dmm=smin(ii)-s(ii)
               CALL vkt_limiter(dels,dpp,dmm,eps2,psi)
               spsi(ii)=DMIN1(psi,spsi(ii))
!
         IF(kk.le.ncell_fluid) THEN
               dels= dsdx(kk,1)*dxfc_non_k(i,1) &
                    +dsdx(kk,2)*dxfc_non_k(i,2)
               dpp=smax(kk)-s(kk)
               dmm=smin(kk)-s(kk)
               CALL vkt_limiter(dels,dpp,dmm,eps2,psi)
               spsi(kk)=DMIN1(psi,spsi(kk))
         ENDIF
            ENDDO
!
!........The rest
            DO nf_number=2,8
               istart=istart_nf(1,nf_number)
               len   =istart_nf(2,nf_number)
               DO i=1,len  
                  i1=istart+i
                  ii=left_nf(i1)
!
               delh=vol(ii)**(1.d0/2.d0)
               IF(vkt.gt.0.0d0)THEN
                  tmp=vkt*delh
                  eps2=tmp*tmp*tmp
               ELSE
                 eps2=0.0d0
               ENDIF
                  dels= dsdx(ii,1)*dxfc_nf(i1,1) &
                       +dsdx(ii,2)*dxfc_nf(i1,2)
                  dpp=smax(ii)-s(ii)
                  dmm=smin(ii)-s(ii)
                  CALL vkt_limiter(dels,dpp,dmm,eps2,psi)
                  spsi(ii)=DMIN1(psi,spsi(ii))
               ENDDO
            ENDDO
         ELSE
            nf_number=0
            len   =istart_nf(2,nf_number)
            DO i=1,len  
               i1=istart+i
               ii=left_nf(i1)
               kk=right_non(i)
!
         delh=vol(ii)**(1.d0/3.d0)
         IF(vkt.gt.0.0d0)THEN
            tmp=vkt*delh
            eps2=tmp*tmp*tmp
         ELSE
            eps2=0.0d0
         ENDIF
!
               dels= dsdx(ii,1)*dxfc_nf(i1,1) &
                    +dsdx(ii,2)*dxfc_nf(i1,2) &
                    +dsdx(ii,3)*dxfc_nf(i1,3)
               dpp=smax(ii)-s(ii)
               dmm=smin(ii)-s(ii)
               CALL vkt_limiter(dels,dpp,dmm,eps2,psi)
               spsi(ii)=DMIN1(psi,spsi(ii))
!
         IF(kk.le.ncell_fluid) THEN
               dels= dsdx(kk,1)*dxfc_non_k(i,1) &
                    +dsdx(kk,2)*dxfc_non_k(i,2) &
                    +dsdx(kk,3)*dxfc_non_k(i,3)
               dpp=smax(kk)-s(kk)
               dmm=smin(kk)-s(kk)
               CALL vkt_limiter(dels,dpp,dmm,eps2,psi)
               spsi(kk)=DMIN1(psi,spsi(kk))
         ENDIF
            ENDDO
!
!........The rest
            DO nf_number=2,8
               istart=istart_nf(1,nf_number)
               len   =istart_nf(2,nf_number)
               DO i=1,len  
                  i1=istart+i
                  ii=left_nf(i1)
                  dels= dsdx(ii,1)*dxfc_nf(i1,1) &
                       +dsdx(ii,2)*dxfc_nf(i1,2) &
                       +dsdx(ii,3)*dxfc_nf(i1,3)
                  dpp=smax(ii)-s(ii)
                  dmm=smin(ii)-s(ii)
                  CALL vkt_limiter(dels,dpp,dmm,eps2,psi)
                  spsi(ii)=DMIN1(psi,spsi(ii))
               ENDDO
            ENDDO
         ENDIF 
!
      ENDIF
!
!.....Remove boundary node gradient
!
      DO nf_number=2,8
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
         DO i=1,len  
            i1=istart+i
            ii=left_nf(i1)
!            spsi(ii)=0.0d0
         ENDDO
      ENDDO
!
!.....Finally, apply the limiter
!
      IF(ndim.eq.2) THEN
         DO i=1,ncell_fluid
            fs=faclim*spsi(i)
            dsdx(i,1)=fs*dsdx(i,1)
            dsdx(i,2)=fs*dsdx(i,2)
         ENDDO
      ELSE
         DO i=1,ncell_fluid
            fs=faclim*spsi(i)
            dsdx(i,1)=fs*dsdx(i,1)
            dsdx(i,2)=fs*dsdx(i,2)
            dsdx(i,3)=fs*dsdx(i,3)
         ENDDO
      ENDIF
!
      RETURN
      END SUBROUTINE grad_limiter
!DEC$ ATTRIBUTES INLINE :: barth_limiter
!
      SUBROUTINE barth_limiter(delq,dpp,dmm,phi)
!
      IMPLICIT NONE
!
      REAL(8) :: delq,dpp,dmm,phi
!
      IF(DABS(delq).lt.1.d-12) delq=0.d0
      IF(DABS(dpp).lt.1.d-12) dpp=0.d0
      IF(DABS(dmm).lt.1.d-12) dmm=0.d0

      IF(delq.gt.0.d0) THEN
         phi=DMIN1(1.d0,dpp/delq)
      ELSEIF(delq.lt.0.d0) THEN
         phi=DMIN1(1.d0,dmm/delq)
      ELSE
         phi=1.d0
      ENDIF
!
      RETURN
      END SUBROUTINE barth_limiter
!DEC$ ATTRIBUTES INLINE :: vkt_limiter
!
      SUBROUTINE vkt_limiter(delq,dpp,dmm,eps2,phi)
!
      IMPLICIT NONE
!
      REAL(8) :: delq,dpp,dmm,eps2,phi
      REAL(8) :: ddd,dd2,dq2,dqd
!
      IF(DABS(delq).lt.1.d-12) delq=0.d0
      IF(DABS(dpp).lt.1.d-12) dpp=0.d0
      IF(DABS(dmm).lt.1.d-12) dmm=0.d0
!
      IF(delq.gt.1.d-12) THEN
         ddd=dpp
         dd2=ddd*ddd ; dq2=delq*delq
         dqd=ddd*delq
         phi=(dd2+2.d0*dqd+eps2)/(dd2+dqd+2.d0*dq2+eps2)
      ELSEIF(delq.lt.-1.d-12) THEN
         ddd=dmm
         dd2=ddd*ddd ; dq2=delq*delq
         dqd=ddd*delq
         phi=(dd2+2.d0*dqd+eps2)/(dd2+dqd+2.d0*dq2+eps2)
      ELSE
         phi=1.d0
      ENDIF
!
      RETURN
      END SUBROUTINE vkt_limiter
