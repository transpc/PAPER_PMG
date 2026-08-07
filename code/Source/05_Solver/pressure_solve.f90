!
      SUBROUTINE pressure_solve(ic,ip,pcnv)
!
!     This routine establish pressure matrix which is solved using conjugate gradient solvers
!     For parallel computing, CZR format is applied for pressure matrix.
!
      USE Zinterface
      USE Zmpi          , ONLY: ncell_fp,au,ia_a,ju_a,maxmt
      USE Zzone         , ONLY: ncell_fluid
      USE Zcore         , ONLY: np,myrank_cupid => myrank
      USE Zparam        , ONLY: ndim
      USE Zvec_param    , ONLY: nf_nonk,nf_non
      USE Znum_cell     , ONLY: istart_nf, &
                                nf_number_nb,lens,                           &
                                right_nb_k,istart_nfs,nf_number_id,istart_nf
      USE Zvec_index    , ONLY: left_nf,right_non
      USE Zconst1       , ONLY: parallel,fric_face,nd_face
      USE Zconst2       , ONLY: dt
      USE Zgradoption   , ONLY: iter_grad,non_orth,non_orth_iter,grav_grad
      USE Zbc_index     , ONLY: npb
      USE Zbicg         , ONLY: eps_bicg,max_bicg,min_bicg,relax_p
      USE Zb_condition  , ONLY: pbnd
      USE Zcoord3       , ONLY: volpr
      USE Zmars         , ONLY: n_marsbc,i3invtbl_tmp,ppcup_tmp
      USE Zpress_coeff  , ONLY: coefp_l,coefp_g,coefp_d,coefm_l,coefm_g
      USE Ztimecon      , ONLY: smac,dp_max
      USE Zuserdefined  , ONLY: udfl_outlet_press_user,MG_solver
      USE Zbc_index     , ONLY: i_horizontal_outlet
      USE c3com_cupid   , ONLY: i3invtbl
      USE Zpress        , ONLY: p,pp,dpdx,flag
      USE Zvector       , ONLY: vl_n,vg_n,vd_n,vl_o,vg_o,vd_o,vl_f_non,vg_f_non
      USE Zvec_major    , ONLY: flux_l_nf,flux_g_nf,flux_d_nf
      USE Zvec_geo      , ONLY: xn_nf,sv_nf,sap_nf,dji_x_nf,djia_nf, &
                                xloc_m_non_i,xloc_m_non_k,           &
                                f0,f1,fac1_non,fac_non,              &
                                fac1_non_k,sv_non_k
      USE Zmcp         
      USE Zrv_model     , ONLY: rv_choke,rv_mcp,rv_valve
      USE MD_matrix     , ONLY: u
      USE MD_MG_INDEX, ONLY: icase_mg
!
      IMPLICIT NONE
!      
!DEC$IF defined (mpi_flag)
      INCLUDE 'mpif.h'
!DEC$ENDIF
      INCLUDE '../10_LinkToMARS/c3com.h' 
!
!.....Input
      INTEGER :: ic,ip
!.....Output
      LOGICAL :: pcnv
!.....Local variables
      INTEGER :: i,j,k
      INTEGER :: ii,kk,ix,iter
      INTEGER :: nv,nf_number,istart0,istart,len,i0,i1
      INTEGER :: idx,mm
      INTEGER :: izone=0
      REAL(8) :: dp,cf,tmp,eps
      REAL(8) :: svolli,svolgi,svoldi
      REAL(8) :: epsFactor,eps2,dpx
      REAL(8) :: grdPx1,grdPx2,grdPx3
      REAL(8) :: dpi,dpj
      REAL(8) :: vl_fr,vg_fr
      REAL(8) :: a_l,a_g,a_d
      REAL(8) :: t1,t2,t_tmp,t3,t4
      REAL(8), SAVE :: t_sum=0.d0
      REAL(8), SAVE :: t_sum1=0.d0
!.....Local arrays
      REAL(8),DIMENSION(ncell_fluid,ndim) :: vg1_n,vl1_n
      REAL(8),DIMENSION(ncell_fluid) :: poiss_diag,src
!.....Local vector arrays
      REAL(8),DIMENSION(nf_nonk+nf_non,ndim) :: fpl_non,fpg_non
      REAL(8),DIMENSION(nf_non) :: poiss_non_i
      REAL(8),DIMENSION(nf_nonk) :: poiss_non_k
      REAL(8),DIMENSION(nf_non) :: vl_fr_non,vg_fr_non
!.....Local allocatable arrays
      REAL(8),DIMENSION(:),ALLOCATABLE :: ppp
!
!.....Initialize minimum residual
!
      epsFactor=1.0d0
!
      IF(ic.eq.1.and.ip.eq.0)THEN
         eps=1.0d-3
      ELSE
         eps=eps_bicg
      ENDIF
!
!.....Set the boundary pressure
!
      IF(udfl_outlet_press_user.eq..false. .and. i_horizontal_outlet.le.0)THEN
         DO i=1,ncell_fluid
            IF(npb(i).gt.0) p(i)=pbnd(npb(i))
         ENDDO
      ENDIF
!
!.....Initialization for pararell computing
!
      IF(parallel.ge.1)THEN
!
!........Communicate pressure and pressure gradient for matrix setup
!
         IF(np.gt.1) THEN
            CALL communicate_1d(coefp_g, &
                                coefp_l, &
                                coefp_d)
            CALL communicate_2d(vg_n, &
                                vl_n, &
                                vd_n)
         ENDIF
!
      ENDIF 
!
!.....Set pressure matrix: ic=0 for SMAC series
!
      IF(ic.eq.0)THEN
         CALL pressure_matrix_smac(poiss_diag,poiss_non_i,poiss_non_k,src)
      ELSE
         CALL pressure_matrix(poiss_diag,poiss_non_i,poiss_non_k,src)
      ENDIF      
!
!.....No solve for npb > 0
!
      DO i=1,ncell_fluid
         IF(npb(i).ne.0)THEN
            poiss_diag(i)=1.0d0
            src(i)=0.0d0
         ENDIF
      ENDDO
!
!.....Build directly solver array CSR format
!
      CALL csr_build_a(poiss_diag,poiss_non_i,poiss_non_k)
!
!.....No solve for npb > 0
!
      DO i=1,ncell_fluid
         IF(npb(i).gt.0)THEN
            DO j=ia_a(i),ju_a(i)-1
               au(j)=0.d0
            ENDDO
            DO j=ju_a(i)+1,ia_a(i+1)-1
               au(j)=0.d0
            ENDDO
         ENDIF
      ENDDO
!
!.....Select parallel matrix solver
!
!DEC$IF defined (mpi_flag)
      t1 = MPI_Wtime()
!DEC$ENDIF
      IF(parallel.eq.0)THEN
!
!........BiCG solver for serial calculation
!
         flag=1
         IF(smac.eq.3) eps=1.0d-10
         eps2=eps*epsFactor
         pp=0.d0
         CALL bicgstab(eps2,max_bicg,min_bicg,pp,src,flag)
!
      ELSEIF(parallel.eq.1)THEN
!
! copy data to MG solver! (matrix and source)
         IF (MG_solver) THEN
!           write(*,*)'mg-1'
            u=0.d0
! test-MG
!!DEC$IF defined (mpi_flag)
!      t3 = MPI_Wtime()
!!DEC$ENDIF

            CALL assemble_FVM(1,maxmt, src, au, poiss_diag)!,poiss_csr)

!!DEC$IF defined (mpi_flag)
!      t4 = MPI_Wtime()
       
!         IF (myrank_cupid .eq. 0) WRITE (101, *) 'assmem', t4-t3
   
!!DEC$ENDIF


            CALL SOLVE_GMG(1)
            pp = u
         ELSE
         CALL cupid_solvers(epsFactor,poiss_diag,non_orth,src,pp,izone)
         ENDIF
      ELSE
         STOP 'Check Matrix Solver!'
      ENDIF  
!DEC$IF defined (mpi_flag)
      t2 = MPI_Wtime()
      t_tmp = t_sum+t2-t1
      CALL allreduce_max_r1(t_tmp, t_sum)
      IF (MG_solver .eq. .false.) THEN
          IF (myrank_cupid .eq. 0) WRITE (200+myrank_cupid, *), t_sum
      ELSE
         IF (myrank_cupid .eq. 0) WRITE (200+200+myrank_cupid, *), t_sum
      END IF
!DEC$ENDIF
!
!.....Check maximum pressure correction and convergence
!
      dp_max=0.0d0
      DO i=1,ncell_fluid
         dp_max=max(dp_max,abs(pp(i)))
      ENDDO
      IF(np.gt.1) CALL allreducei_max_r1(dp_max)
!
      IF(relax_p.gt.0.0d0) THEN
         DO i=1,ncell_fp
            pp(i)=(1.0d0-relax_p)*pp(i)
         ENDDO
      ENDIF
! 
!.....Calculate pressure correction gradient
! 
      IF(iter_grad.gt.1) THEN
         CALL grad_pressK2(pp,dpdx,1)
      ELSE
         IF(grav_grad.ge.2)THEN
            CALL grad_frink(pp,dpdx,2)
         ELSE
            CALL grad_press(pp,dpdx,1)
         ENDIF
      ENDIF
      IF(np.gt.1) CALL communicate_2d(dpdx)
!
!.....Pressure corection for non-orthogonal grid
!
      ! test-MG
!        write(*,*)'non-or',non_orth,non_orth_iter

      IF(non_orth.gt.0)THEN
         ALLOCATE(ppp(ncell_fp))
         DO iter=1,non_orth_iter
            CALL non_orthogonal_src(src,poiss_non_i,poiss_non_k)
            DO i=1,ncell_fluid
               IF(npb(i).gt.0) src(i)=0.0d0
            ENDDO
            IF(parallel.EQ.0) THEN
               flag=1
               eps2=eps*epsFactor 
               CALL bicgstab(eps_bicg,max_bicg,min_bicg,ppp,src,flag)
            ELSEIF(parallel.eq.1) THEN

      !DEC$IF defined (mpi_flag)
         t3 = MPI_Wtime()
         !DEC$ENDIF

               IF (MG_solver) THEN
!                write(*,*)'MG-2'
                  CALL assemble_FVM(icase_mg,maxmt,src,au,poiss_diag)!,poiss_csr)
                  CALL SOLVE_GMG(icase_mg)
                  ppp = u
               ELSE
               CALL cupid_solve_non_orth(epsFactor,poiss_diag,src,ppp,izone)
               ENDIF

      !DEC$IF defined (mpi_flag)
         t4 = MPI_Wtime()
         !DEC$ENDIF

            ENDIF
            IF(relax_p.gt.0.0d0) THEN
               DO i=1,ncell_fp
                  ppp(i)=pp(i)+(1.0d0-relax_p)*ppp(i)
               ENDDO
            ELSE
               DO i=1,ncell_fp
                  ppp(i)=pp(i)+ppp(i)
               ENDDO
            ENDIF
            IF(grav_grad.ge.2)THEN
               CALL grad_frink(ppp,dpdx,2)
            ELSE
               CALL grad_pressK2(ppp,dpdx,1)
            ENDIF
            IF(np.gt.1) CALL communicate_2d(dpdx)
         ENDDO
         DO i=1,ncell_fp
            pp(i)=ppp(i)
         ENDDO
         DEALLOCATE(ppp)
      ENDIF

      !DEC$IF defined (mpi_flag)
       t_tmp = t_sum1 +t4-t3+t2-t1

         CALL allreduce_max_r1(t_tmp, t_sum1)
         IF (MG_solver .eq. .false.) THEN
            IF (myrank_cupid .eq. 0) WRITE (201+myrank_cupid, *), t_sum1
         ELSE
            IF (myrank_cupid .eq. 0) WRITE (200+201+myrank_cupid, *),t_sum1
         END IF
         !DEC$ENDIF


!
!.....Initialize variables for MARS code coupling
!
      IF(n_marsbc.gt.0) CALL pickup_pp_cupvols
!
!.....Cells non
      nf_number=0
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      IF(non_orth.eq.0)THEN
         DO i=1,len  
!            IF(npb(left_non(i)).ne.0) CYCLE !Because pressure cell can
!            be left_non(i), where i is face.
           i1=istart+i
           ii=left_nf(i1)
           kk=right_non(i)
           cf=dt*sap_nf(i1)
!                    
           dp=pp(kk)-pp(ii)
           svolli=f1(i)*coefp_l(ii)+f0(i)*coefp_l(kk)
           svolgi=f1(i)*coefp_g(ii)+f0(i)*coefp_g(kk)
           svoldi=f1(i)*coefp_d(ii)+f0(i)*coefp_d(kk)
           flux_l_nf(i1)=flux_l_nf(i1)-cf*svolli*dp
           flux_g_nf(i1)=flux_g_nf(i1)-cf*svolgi*dp
           flux_d_nf(i1)=flux_d_nf(i1)-cf*svoldi*dp
         ENDDO
!          
      ELSEIF(non_orth.eq.1)THEN
         IF(ndim.eq.2)THEN
             DO i=1,len  
!                IF(npb(left_non(i)).ne.0) CYCLE
                i1=istart+i
                ii=left_nf(i1)
                kk=right_non(i)
                cf=dt*sap_nf(i1)
!                      
                dp=pp(kk)-pp(ii)
                grdPx1=0.5d0*(dpdx(ii,1)+dpdx(kk,1))
                grdPx2=0.5d0*(dpdx(ii,2)+dpdx(kk,2))
                dpx=grdPx1*xn_nf(i1,1)+grdPx2*xn_nf(i1,2)
                dp=dp+dpx*djia_nf(i1)
                dpx=grdPx1*dji_x_nf(i1,1)+grdPx2*dji_x_nf(i1,2)
                dp=dp-dpx
                svolli=fac1_non(i)*coefp_l(ii)+fac_non(i)*coefp_l(kk)
                svolgi=fac1_non(i)*coefp_g(ii)+fac_non(i)*coefp_g(kk)
                svoldi=fac1_non(i)*coefp_d(ii)+fac_non(i)*coefp_d(kk)
                flux_l_nf(i1)=flux_l_nf(i1)-cf*svolli*dp
                flux_g_nf(i1)=flux_g_nf(i1)-cf*svolgi*dp
                flux_d_nf(i1)=flux_d_nf(i1)-cf*svoldi*dp
             ENDDO
         ELSE
             DO i=1,len  
!                IF(npb(left_non(i)).ne.0) CYCLE
                i1=istart+i
                ii=left_nf(i1)
                kk=right_non(i)
                cf=dt*sap_nf(i1)
!                         
                dp=pp(kk)-pp(ii)
                grdPx1=0.5d0*(dpdx(ii,1)+dpdx(kk,1))
                grdPx2=0.5d0*(dpdx(ii,2)+dpdx(kk,2))
                grdPx3=0.5d0*(dpdx(ii,3)+dpdx(kk,3))
                dpx=grdPx1*xn_nf(i1,1)+grdPx2*xn_nf(i1,2)+grdPx3*xn_nf(i1,3)
                dp=dp+dpx*djia_nf(i1)
                dpx=grdPx1*dji_x_nf(i1,1)+grdPx2*dji_x_nf(i1,2)+grdPx3*dji_x_nf(i1,3)
                dp=dp-dpx
                svolli=fac1_non(i)*coefp_l(ii)+fac_non(i)*coefp_l(kk)
                svolgi=fac1_non(i)*coefp_g(ii)+fac_non(i)*coefp_g(kk)
                svoldi=fac1_non(i)*coefp_d(ii)+fac_non(i)*coefp_d(kk)
                flux_l_nf(i1)=flux_l_nf(i1)-cf*svolli*dp
                flux_g_nf(i1)=flux_g_nf(i1)-cf*svolgi*dp
                flux_d_nf(i1)=flux_d_nf(i1)-cf*svoldi*dp
             ENDDO
         ENDIF
!          
      ELSEIF(non_orth.eq.2)THEN
         IF(ndim.eq.2)THEN
             DO i=1,len  
!                IF(npb(left_non(i)).ne.0) CYCLE
                i1=istart+i
                ii=left_nf(i1)
                kk=right_non(i)
                cf=dt*sap_nf(i1)
!                                
                dp=pp(kk)-pp(ii)
                dpi=dpdx(ii,1)*xloc_m_non_i(i,1)+dpdx(ii,2)*xloc_m_non_i(i,2)
                dpj=dpdx(kk,1)*xloc_m_non_k(i,1)+dpdx(kk,2)*xloc_m_non_k(i,2)
                dp=dp+dpj-dpi
                svolli=fac1_non(i)*coefp_l(ii)+fac_non(i)*coefp_l(kk)
                svolgi=fac1_non(i)*coefp_g(ii)+fac_non(i)*coefp_g(kk)
                svoldi=fac1_non(i)*coefp_d(ii)+fac_non(i)*coefp_d(kk)
                flux_l_nf(i1)=flux_l_nf(i1)-cf*svolli*dp
                flux_g_nf(i1)=flux_g_nf(i1)-cf*svolgi*dp
                flux_d_nf(i1)=flux_d_nf(i1)-cf*svoldi*dp
             ENDDO
         ELSE
             DO i=1,len  
!                IF(npb(left_non(i)).ne.0) CYCLE
                i1=istart+i
                ii=left_nf(i1)
                kk=right_non(i)
                cf=dt*sap_nf(i1)
!                        
                dp=pp(kk)-pp(ii)
                dpi=dpdx(ii,1)*xloc_m_non_i(i,1)+dpdx(ii,2)*xloc_m_non_i(i,2)+dpdx(ii,3)*xloc_m_non_i(i,3)
                dpj=dpdx(kk,1)*xloc_m_non_k(i,1)+dpdx(kk,2)*xloc_m_non_k(i,2)+dpdx(kk,3)*xloc_m_non_k(i,3)
                dp=dp+(dpj-dpi)
                svolli=fac1_non(i)*coefp_l(ii)+fac_non(i)*coefp_l(kk)
                svolgi=fac1_non(i)*coefp_g(ii)+fac_non(i)*coefp_g(kk)
                svoldi=fac1_non(i)*coefp_d(ii)+fac_non(i)*coefp_d(kk)
                flux_l_nf(i1)=flux_l_nf(i1)-cf*svolli*dp
                flux_g_nf(i1)=flux_g_nf(i1)-cf*svolgi*dp
                flux_d_nf(i1)=flux_d_nf(i1)-cf*svoldi*dp
             ENDDO
         ENDIF
      ENDIF 
!
!.....valve model
!      
      IF(rv_valve.eq.1) CALL valve_model_pressure_solve(non_orth)
!     
!.....MARS interface
!
      nf_number=1
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
      IF(ic.eq.0)THEN
!DIR$ SIMD
         DO i=1,len  
            i1=istart+i
            ii=left_nf(i1)
            idx=i3invtbl(i)
            cf=dt*c3area(1,idx)
            tmp=0.0d0
            DO mm=1,n_marsbc
               k=i3invtbl_tmp(mm)
               tmp=tmp+c3yeta(1,idx,k)*ppcup_tmp(k) !cupid_mars_debug
            ENDDO
            flux_l_nf(i1)=flux_l_nf(i1)+cf*c3betaf(1,idx)*(pp(ii)-c3xi(1,idx)-tmp)
            flux_g_nf(i1)=flux_g_nf(i1)+cf*c3betag(1,idx)*(pp(ii)-c3xi(1,idx)-tmp)
            flux_d_nf(i1)=flux_d_nf(i1)+cf*c3betaf(1,idx)*(pp(ii)-c3xi(1,idx)-tmp)
            c3vg(1,idx)=c3alphg(1,idx)+c3betag(1,idx)*(pp(ii)-c3xi(1,idx)-tmp)
            c3vl(1,idx)=c3alphf(1,idx)+c3betaf(1,idx)*(pp(ii)-c3xi(1,idx)-tmp)
         ENDDO
      ELSE
!DIR$ SIMD
         DO i=1,len  
            i1=istart+i
            ii=left_nf(i1)
            idx=i3invtbl(i)
            cf=-dt*sap_nf(i1)
            tmp=0.0d0
            DO mm=1,n_marsbc
               k=i3invtbl_tmp(mm)
               tmp=tmp+c3yeta(1,idx,k)*ppcup_tmp(k) !cupid_mars_debug
            ENDDO
            flux_l_nf(i1)=flux_l_nf(i1)+cf*c3betaf(1,idx)*(pp(ii)-c3xi(1,idx)-tmp)
            flux_g_nf(i1)=flux_g_nf(i1)+cf*c3betag(1,idx)*(pp(ii)-c3xi(1,idx)-tmp)
            flux_d_nf(i1)=flux_d_nf(i1)+cf*c3betaf(1,idx)*(pp(ii)-c3xi(1,idx)-tmp)
         ENDDO
      ENDIF
!
!.....fluxBC: choke model, mcp model, valve model (flux(n+1)=flux(choked), u(*)=u(choked))
! 
      IF(rv_choke.eq.1.or.rv_mcp.eq.1.or.rv_valve.eq.1) CALL fluxBC_flux_update(flux_l_nf,flux_g_nf,flux_d_nf)      
!
!.....Correct pressure and velocities
!      
      IF(ndim.eq.2) THEN
         DO i=1,ncell_fluid
            IF(npb(i).ne.0) cycle
            a_l=coefp_l(i)*dt
            a_g=coefp_g(i)*dt
            a_d=coefp_d(i)*dt
            vl_n(i,1)=vl_n(i,1)-a_l*dpdx(i,1)
            vg_n(i,1)=vg_n(i,1)-a_g*dpdx(i,1)
            vd_n(i,1)=vd_n(i,1)-a_d*dpdx(i,1)
            vl_n(i,2)=vl_n(i,2)-a_l*dpdx(i,2)
            vg_n(i,2)=vg_n(i,2)-a_g*dpdx(i,2)
            vd_n(i,2)=vd_n(i,2)-a_d*dpdx(i,2)
         ENDDO
      ELSE
         DO i=1,ncell_fluid
            IF(npb(i).ne.0) cycle
            a_l=coefp_l(i)*dt
            a_g=coefp_g(i)*dt
            a_d=coefp_d(i)*dt
            vl_n(i,1)=vl_n(i,1)-a_l*dpdx(i,1)
            vg_n(i,1)=vg_n(i,1)-a_g*dpdx(i,1)
            vd_n(i,1)=vd_n(i,1)-a_d*dpdx(i,1)
            vl_n(i,2)=vl_n(i,2)-a_l*dpdx(i,2)
            vg_n(i,2)=vg_n(i,2)-a_g*dpdx(i,2)
            vd_n(i,2)=vd_n(i,2)-a_d*dpdx(i,2)
            vl_n(i,3)=vl_n(i,3)-a_l*dpdx(i,3)
            vg_n(i,3)=vg_n(i,3)-a_g*dpdx(i,3)
            vd_n(i,3)=vd_n(i,3)-a_d*dpdx(i,3)
         ENDDO
      ENDIF
!
!.....fluxBC: choke model, mcp model, valve model (flux(n+1)=flux(choked), u(n+1)=u(choked))      
!      
!      IF(rv_choke.eq.1.or.rv_mcp.eq.1.or.rv_valve.eq.1) CALL fluxBC_flux_update(flux_l_nf,flux_g_nf,flux_d_nf)            
!
      DO i=1,ncell_fluid
         IF(p(i)+pp(i).ge.0.0d0) p(i)=p(i)+pp(i)
      ENDDO
!
!.....Interpolation of cell center velocity from face values
!
      IF(fric_face+nd_face.gt.0)THEN
!
!........Build summation info for non
!
         nf_number_nb=0
         nf_number_id(-1)=-1
         nf_number_id(0)=0
         istart_nfs(-1)=0
         istart_nfs(0)=istart_nfs(-1)+nf_nonk
         lens         =istart_nfs(0) +nf_non
!
         nv=0
         nf_number=nf_number_id(nv)
         istart0=istart_nfs(nv)
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
         IF(ndim.eq.2) THEN
            DO i=1,len
               i0=istart0+i
               i1=istart+i
               vl_fr_non(i)=vl_f_non(i,1)*dji_x_nf(i1,1)+vl_f_non(i,2)*dji_x_nf(i1,2)
               vg_fr_non(i)=vg_f_non(i,1)*dji_x_nf(i1,1)+vg_f_non(i,2)*dji_x_nf(i1,2)
               vl_fr=fac_non(i)*vl_fr_non(i)
               vg_fr=fac_non(i)*vg_fr_non(i)
               fpl_non(i0,1)=sv_nf(i1,1)*vl_fr
               fpg_non(i0,1)=sv_nf(i1,1)*vg_fr
               fpl_non(i0,2)=sv_nf(i1,2)*vl_fr
               fpg_non(i0,2)=sv_nf(i1,2)*vg_fr
            ENDDO
         ELSE
            DO i=1,len
               i0=istart0+i
               i1=istart+i
               vl_fr_non(i)=vl_f_non(i,1)*dji_x_nf(i1,1)+vl_f_non(i,2)*dji_x_nf(i1,2)+vl_f_non(i,3)*dji_x_nf(i1,3)
               vg_fr_non(i)=vg_f_non(i,1)*dji_x_nf(i1,1)+vg_f_non(i,2)*dji_x_nf(i1,2)+vg_f_non(i,3)*dji_x_nf(i1,3)
               vl_fr=fac_non(i)*vl_fr_non(i)
               vg_fr=fac_non(i)*vg_fr_non(i)
               fpl_non(i0,1)=sv_nf(i1,1)*vl_fr
               fpg_non(i0,1)=sv_nf(i1,1)*vg_fr
               fpl_non(i0,2)=sv_nf(i1,2)*vl_fr
               fpg_non(i0,2)=sv_nf(i1,2)*vg_fr
               fpl_non(i0,3)=sv_nf(i1,3)*vl_fr
               fpg_non(i0,3)=sv_nf(i1,3)*vg_fr
            ENDDO
         ENDIF
!
         nv=-1
         nf_number=nf_number_id(nv)
         len   =istart_nf(2,nf_number)
         IF(ndim.eq.2) THEN
            DO i=1,len
               k=right_nb_k(i)
               vl_fr=-fac1_non_k(i)*vl_fr_non(k)
               vg_fr=-fac1_non_k(i)*vg_fr_non(k)
               fpl_non(i,1)=-sv_non_k(i,1)*vl_fr
               fpg_non(i,1)=-sv_non_k(i,1)*vg_fr
               fpl_non(i,2)=-sv_non_k(i,2)*vl_fr
               fpg_non(i,2)=-sv_non_k(i,2)*vg_fr
            ENDDO
         ELSE
            DO i=1,len
               k=right_nb_k(i)
               vl_fr=-fac1_non_k(i)*vl_fr_non(k)
               vg_fr=-fac1_non_k(i)*vg_fr_non(k)
               fpl_non(i,1)=-sv_non_k(i,1)*vl_fr
               fpg_non(i,1)=-sv_non_k(i,1)*vg_fr
               fpl_non(i,2)=-sv_non_k(i,2)*vl_fr
               fpg_non(i,2)=-sv_non_k(i,2)*vg_fr
               fpl_non(i,3)=-sv_non_k(i,3)*vl_fr
               fpg_non(i,3)=-sv_non_k(i,3)*vg_fr
            ENDDO
         ENDIF
!
         CALL sum_nf_ndim(0,0,ncell_fluid, &
                          fpg_non,vg1_n,   &
                          fpl_non,vl1_n)
!
         IF(ndim.eq.2) THEN
            DO i=1,ncell_fluid
               a_l=coefm_l(i)*volpr(i)
               a_g=coefm_g(i)*volpr(i)
               vl_n(i,1)=vl_n(i,1)+a_l*vl1_n(i,1)
               vl_n(i,2)=vl_n(i,2)+a_l*vl1_n(i,2)
               vg_n(i,1)=vg_n(i,1)+a_g*vg1_n(i,1)
               vg_n(i,2)=vg_n(i,2)+a_g*vg1_n(i,2)
            ENDDO
         ELSE
            DO i=1,ncell_fluid
               a_l=coefm_l(i)*volpr(i)
               a_g=coefm_g(i)*volpr(i)
               vl_n(i,1)=vl_n(i,1)+a_l*vl1_n(i,1)
               vl_n(i,2)=vl_n(i,2)+a_l*vl1_n(i,2)
               vl_n(i,3)=vl_n(i,3)+a_l*vl1_n(i,3)
               vg_n(i,1)=vg_n(i,1)+a_g*vg1_n(i,1)
               vg_n(i,2)=vg_n(i,2)+a_g*vg1_n(i,2)
               vg_n(i,3)=vg_n(i,3)+a_g*vg1_n(i,3)
            ENDDO
         ENDIF
      ENDIF
!
      IF(.not.pcnv) THEN
         DO ix=1,ndim
            DO i=1,ncell_fluid
               vl_n(i,ix)=vl_o(i,ix)
               vg_n(i,ix)=vg_o(i,ix)
               vd_n(i,ix)=vd_o(i,ix)
            ENDDO
         ENDDO
         IF(np.gt.1) CALL communicate_1d(p)
      ENDIF
!
      END SUBROUTINE pressure_solve
