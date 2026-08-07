!
      SUBROUTINE rv_hts_1d
!
      USE Zrv_hts_1d,   ONLY:ncell_hts_1d,nr_1d,nmat_1d,imp_cond_1d,vl_1d,vr_1d,sl_1d,    &
                              sr_1d,bcl_1d,bcr_1d,slw_1d,srw_1d,t_hts_1d,                  &
                              hll_1d,hstl_1d,hspl_1d,hgl_1d,tll_1d,tstl_1d,tspl_1d,tgl_1d, &
                              hlr_1d,hstr_1d,hspr_1d,hgr_1d,tlr_1d,tstr_1d,tspr_1d,tgr_1d, &
                              hfluxl_1d,hfluxr_1d,twl_1d,twr_1d,ig_hts_1d
      USE Zconst2,      ONLY:dt
!
      IMPLICIT NONE
!     local variables
      INTEGER i,k,m
      REAL(8) :: g
!     local arrays
      REAL(8) :: rcp(ncell_hts_1d,nr_1d),cond(ncell_hts_1d,nr_1d)
      REAL(8) :: a(ncell_hts_1d,nr_1d),b(ncell_hts_1d,nr_1d),c(ncell_hts_1d,nr_1d),d(ncell_hts_1d,nr_1d)
!
!
      DO i=1,nr_1d-1
         DO k=1,ncell_hts_1d
            m=ig_hts_1d(k)
!           t_i=0.5d0*(t_hts_1d(k,i)+t_hts_1d(k,i+1))
!           CALL mat_prop(nmat_1d(m),t_i,rcp(k,i),cond(k,i),iOKr,iOKk)
         ENDDO
      ENDDO
      CALL mat_prop_2d(nmat_1d,ig_hts_1d,t_hts_1d,rcp,cond)
!
!........Left wall boundary conditions
!
         i=1
!
      IF(imp_cond_1d.eq.1)THEN      ! Fully implicit
         DO k=1,ncell_hts_1d
            m=ig_hts_1d(k)
            g=vr_1d(m,i)*rcp(k,i)
            a(k,i)=0.0d0
            c(k,i)=-cond(k,i)*sr_1d(m,i)*dt
            d(k,i)=g*t_hts_1d(k,i)
            b(k,i)=g-c(k,i)
         ENDDO !k
      ELSE                       ! Crank-Nicholson
         DO k=1,ncell_hts_1d
            m=ig_hts_1d(k)
            g=vr_1d(m,i)*rcp(k,i)
            a(k,i)=0.0d0
            c(k,i)=-0.5d0*cond(k,i)*sr_1d(m,i)*dt
            d(k,i)=(g+c(k,i))*t_hts_1d(k,i)-c(k,i)*t_hts_1d(k,i+1)
            b(k,i)=g-c(k,i)
         ENDDO !k
      ENDIF
!
      DO k=1,ncell_hts_1d
         m=ig_hts_1d(k)
         IF(bcl_1d(m).eq.1)THEN      ! h + t
            b(k,i)=b(k,i)+(hll_1d(k)+hstl_1d(k)+hspl_1d(k)+hgl_1d(k))*slw_1d(m)*dt
            d(k,i)=d(k,i)+(hll_1d(k)*tll_1d(k)+hstl_1d(k)*tstl_1d(k)+hspl_1d(k)*tspl_1d(k)+hgl_1d(k)*tgl_1d(k))*slw_1d(m)*dt
         ELSEIF(bcl_1d(m).eq.2)THEN ! Constant heat flux
            d(k,i)=d(k,i)+hfluxl_1d(m)*slw_1d(m)*dt
         ELSEIF(bcl_1d(m).eq.3)THEN ! Constant T_wall
            b(k,i)=1.0d0
            c(k,i)=0.0d0
            d(k,i)=twl_1d(m)
         ELSEIF(bcl_1d(m).eq.4)THEN ! Adiabatic wall
         ELSE
            STOP '### Left wall BC is not specified.'
         ENDIF
      ENDDO !k
!
      IF(imp_cond_1d.eq.1)THEN      ! Fully implicit
         DO i=2,nr_1d-1
            DO k=1,ncell_hts_1d
               m=ig_hts_1d(k)
               g=vl_1d(m,i)*rcp(k,i-1)+vr_1d(m,i)*rcp(k,i)
               a(k,i)=-cond(k,i-1)*sl_1d(m,i)*dt
               c(k,i)=-cond(k,i)*sr_1d(m,i)*dt
               b(k,i)=g-a(k,i)-c(k,i)
               d(k,i)=g*t_hts_1d(k,i)
            ENDDO 
         ENDDO
      ELSE                       ! Crank-Nicholson
         DO i=2,nr_1d-1
            DO k=1,ncell_hts_1d
               m=ig_hts_1d(k)
               g=vl_1d(m,i)*rcp(k,i-1)+vr_1d(m,i)*rcp(k,i)
               a(k,i)=-0.5d0*cond(k,i-1)*sl_1d(m,i)*dt
               c(k,i)=-0.5d0*cond(k,i)*sr_1d(m,i)*dt
               b(k,i)=g-a(k,i)-c(k,i)
               d(k,i)=-a(k,i)*t_hts_1d(k,i-1)+(g+a(k,i)+c(k,i))*t_hts_1d(k,i)-c(k,i)*t_hts_1d(k,i+1)
            ENDDO
         ENDDO
      ENDIF
!
!........Right wall boundary conditions
!
         i=nr_1d
!
      IF(imp_cond_1d.eq.1)THEN      ! Fully implicit
         DO k=1,ncell_hts_1d
            m=ig_hts_1d(k)
            g=vl_1d(m,i)*rcp(k,i-1)
            c(k,i)=0.0d0
            a(k,i)=-cond(k,i-1)*sl_1d(m,i)*dt
            d(k,i)=g*t_hts_1d(k,i)
            b(k,i)=g-a(k,i)
         ENDDO 
      ELSE                       ! Crank-Nicholson
         DO k=1,ncell_hts_1d
            m=ig_hts_1d(k)
            g=vl_1d(m,i)*rcp(k,i-1)
            c(k,i)=0.0d0
            a(k,i)=-0.5d0*cond(k,i-1)*sl_1d(m,i)*dt
            d(k,i)=-a(k,i)*t_hts_1d(k,i-1)+(g+a(k,i))*t_hts_1d(k,i)
            b(k,i)=g-a(k,i)
         ENDDO 
      ENDIF
!
      DO k=1,ncell_hts_1d
         m=ig_hts_1d(k)
         IF(bcr_1d(m).eq.1)THEN      ! h + t
            b(k,i)=b(k,i)+(hlr_1d(k)+hstr_1d(k)+hspr_1d(k)+hgr_1d(k))*srw_1d(m)*dt
            d(k,i)=d(k,i)+(hlr_1d(k)*tlr_1d(k)+hstr_1d(k)*tstr_1d(k)+hspr_1d(k)*tspr_1d(k)+hgr_1d(k)*tgr_1d(k))*srw_1d(m)*dt         
         ELSEIF(bcr_1d(m).eq.2)THEN ! Constant heat flux
            d(k,i)=d(k,i)+hfluxr_1d(m)*srw_1d(m)*dt         
         ELSEIF(bcr_1d(m).eq.3)THEN ! Constant T_wall
            a(k,i)=0.0d0
            b(k,i)=1.0d0
            d(k,i)=twr_1d(m)
         ELSEIF(bcr_1d(m).eq.4)THEN ! Adiabatic wall
         ELSE
            STOP '### Right wall BC is not specified.'
         ENDIF
      ENDDO 
!
!........Solve tridiagonal matrix
!
         CALL rv_tdiag_2d(a,b,c,d,ncell_hts_1d,nr_1d)
!
      DO i=1,nr_1d
         DO k=1,ncell_hts_1d
            t_hts_1d(k,i)=d(k,i)
         ENDDO
      ENDDO 
!
!
      RETURN
      END SUBROUTINE rv_hts_1d
