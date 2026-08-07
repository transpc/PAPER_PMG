!
      SUBROUTINE udfn_heat_tube
!
      USE Zzone       ,ONLY: ncell_fluid
      USE Zcore       ,ONLY: myrank,np
      USE Ztimecon    ,ONLY: time
      USE Zconst2     ,ONLY: dt
      USE Zsg         ,ONLY: n_group,nr_tube,t_tube,p_tube,vol_tube,dr_tube,dr_tube2,t_1d,q_pri,q_sd, &
                             h_tube,vol_1d,ng_tube,htc_pr,htcl_sec,q_rated,max_1d,q_mult_pri,      &
                             time_sg_heat_tune,t_wall,htcg_sec,htcb_sec,tl_sg,tg_sg,ts_sg,igr,j1d
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER g,m,k,i,iOKr,iOKk
      LOGICAL,SAVE :: initial=.true.
      REAL(8) a1,a2,htc_wall1,tmp,cd2,htc
      REAL(8) q_pri_tot,q_sd_tot,q_pri_tot1,q_sd_tot1
!
!.....Local arrays
      REAL(8) :: temp(2)
      REAL(8),ALLOCATABLE::a(:),b(:),c(:),d(:),cond_tube(:),rcp_tube(:),t_old(:)
      REAL(8),ALLOCATABLE::q_pri_tmp(:,:),q_sd_tmp(:,:)
!
      ALLOCATE(a(nr_tube),b(nr_tube),c(nr_tube),d(nr_tube))
      ALLOCATE(cond_tube(nr_tube),rcp_tube(nr_tube),t_old(nr_tube))
      ALLOCATE(q_pri_tmp(n_group,max_1d),q_sd_tmp(n_group,max_1d))
!
      IF(initial.and.myrank.eq.0)THEN
         OPEN(458,file='ht_area_mult.out')
         initial=.false.
      ENDIF
!
      q_pri_tot=0.0d0
      q_sd_tot=0.0d0
      q_pri_tot1=0.0d0
      q_sd_tot1=0.0d0
!
      DO i=1,ncell_fluid
         g=igr(i)
         IF(g.eq.0) CYCLE
         m=j1d(i)
!
         a(:)=0.0d0
         b(:)=0.0d0
         c(:)=0.0d0
         d(:)=0.0d0
         t_old(:)=t_tube(:,g,m)
!
         DO k=1,nr_tube
            CALL mat_prop(3,t_tube(k,g,m),rcp_tube(k),cond_tube(k),iOKr,iOKk)
         ENDDO
!
  100    CONTINUE 
!
         DO k=2,nr_tube-1
            a(k)=-p_tube(k)*(cond_tube(k-1)+cond_tube(k))/dr_tube2
            c(k)=-p_tube(k+1)*(cond_tube(k+1)+cond_tube(k))/dr_tube2
            a1=vol_tube(k)*rcp_tube(k)/dt
            b(k)=a1-a(k)-c(k)
            d(k)=a1*t_old(k)
         ENDDO
!
!........Primary surface
!
         a(1)=0.0d0
         c(1)=-p_tube(2)*(cond_tube(2)+cond_tube(1))/dr_tube2
         a1=vol_tube(1)*rcp_tube(1)/dt
         a2=2.0d0*cond_tube(1)*htc_pr(i)/(2.0d0*cond_tube(1)+htc_pr(i)*dr_tube)*p_tube(1)
         b(1)=a1-c(1)+a2
         d(1)=a1*t_old(1)+a2*t_1d(g,m)
!
!........Secondary surface
!
         c(nr_tube)=0.0d0
         a(nr_tube)=-p_tube(nr_tube)*(cond_tube(nr_tube-1)+cond_tube(nr_tube))/dr_tube2
         a1=vol_tube(nr_tube)*rcp_tube(nr_tube)/dt
         cd2=2.0d0*cond_tube(nr_tube)
         htc=htcl_sec(i)+htcg_sec(i)+htcb_sec(i)
         tmp=cd2+htc*dr_tube
         a2=p_tube(nr_tube+1)/tmp
         b(nr_tube)=a1-a(nr_tube)+cd2*htc*a2
         d(nr_tube)=a1*t_old(nr_tube)+   &
                    cd2*(htcl_sec(i)*tl_sg(g,m)+htcg_sec(i)*tg_sg(g,m)+htcb_sec(i)*ts_sg(g,m))*a2
!
         CALL tdiag(a,b,c,d,nr_tube)
!
         t_tube(:,g,m)=d(:)
!
!........Set primary surface heat flux
!
         htc_wall1=2.0d0*cond_tube(1)*htc_pr(i)/(2.0d0*cond_tube(1)+htc_pr(i)*dr_tube)
         q_pri_tmp(g,m)=htc_wall1*p_tube(1)*h_tube(g,m)*(t_tube(1,g,m)-t_1d(g,m))*ng_tube(g)
!
!........Set secondary side wall temperature
!
         t_wall(g,m)=cd2*t_tube(nr_tube,g,m)+(htcl_sec(i)*tl_sg(g,m)+htcb_sec(i)*ts_sg(g,m)+htcg_sec(i)*tg_sg(g,m))*dr_tube
         t_wall(g,m)=t_wall(g,m)/tmp
!
!........Calculate totla heat transfer rate
!
         q_pri_tot=q_pri_tot+q_pri_tmp(g,m)
         q_sd_tot=q_sd_tot+q_sd(g,m)
!
      ENDDO
!
      IF(np.gt.1)THEN
         temp(1)=q_pri_tot
         temp(2)=q_sd_tot
         CALL allreducei_r(temp,2)
         q_pri_tot=temp(1)
         q_sd_tot =temp(2)
      ENDIF
!
!.....Calculate heat transfer tunning factor
!
      IF(time.lt.time_sg_heat_tune)THEN
         q_mult_pri=-q_rated/q_pri_tot
         q_mult_pri=MAX(0.1d0,q_mult_pri)
      ENDIF
!
!.....Set the heat sink and source for the primary and secondary coolant
!
      DO i=1,ncell_fluid
         g=igr(i)
         IF(g.eq.0) CYCLE
         m=j1d(i)
         q_pri(i)=q_mult_pri*q_pri_tmp(g,m)/vol_1d(g,m)
      ENDDO
!
      DEALLOCATE(a,b,c,d,cond_tube,rcp_tube,t_old,q_pri_tmp,q_sd_tmp)
!
      IF(myrank.eq.0) WRITE(458,*) time,q_mult_pri
!
      END SUBROUTINE udfn_heat_tube
