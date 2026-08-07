!
!     Those routines provide a user friendly interface for the summation of neighbors involved in a cell
!     The neighbors are arranged in stride 1 for every nf number except for the kk where k<i
!     The block of associated neighbors is accessed with ia_nb in CSR style fashion 
!     For the kk processing each element is extracted with right_nb_k index out of the non vector space
!
!     The calling interface is designed to: 
!          1: make each subroutine independent of the number of vectors processed and shared
!          2: share summation routines
!          3: use contiguous storage no holes no memory waisted
!          4: compute stride 1 for efficiency except for the k<i where vector data is not stored 
!
!             Example: summation process non,inl,fsw vectors for 2 variables
!
!             Set input: x1_input_nf(nf_non+nf_inl+nf_fsw) <==contiguous
!             Set input: x2_input_nf(nf_non+nf_inl+nf_fsw) <==contiguous
!
!             Compute input in vector mode
! 
!             Set control parameters in Znum_cell
!             nf_number_nb=3
!             nf_number_id(1)=0  <==non 
!             nf_number_id(2)=2  <==inl
!             nf_number_id(3)=5  <==fsw
!             istart_nfs(1)=0
!             istart_nfs(2)=istart_nfs(1)+nf_non
!             istart_nfs(3)=istart_nfs(1)+nf_inl
!             lens=istart_nfs(3)+nf_fsw           <==total size compressed vector
!
!             to compute input vector for fsw

!           nv=3
!           nf_number=nf_number_id(nv)
!           istart0=istart_nfs(nv)
!           istart=istart_nf(1,nf_number)
!           len   =istart_nf(2,nf_number)
!           DO i=1,len
!              i0=istart0+i
!              i1=istart+i
!              ii=left_nf(i1)
!              kk=right_non(i)
!              x2_input_nf(i0)=sap_nf(i1)  <===i0 adresses compressed vector 
!                                          <===i1 adresses vector space
!
!             CALL sum_nf(isign,zero,sym,          &
!                         x1_input_nf,x2_input_nf, &
!                         y1_output,y2_output)
!
!--------------------------------------------------------------------------
!     sum_nf(zero,sym,       summation max 8 variables
!                            zero:   0 zero the output first
!                                    1 cumulative output
!                            sym:    0 non symetric k
!                                    1 symetric k, special path for k add
!                                   -1 symetric k, special path for k substarct
!
!     sum_nf_ndim            summation of ndim max 4 variables
!     sum_nf_ndim2           summation of ndim,ndim for 1 variable
!     sum_nf_solid           summation for 1 variable for solid processing
!
!     max_nf                 get the maximum for all for 1 variable 
!     sum_nf0_idc            summation for non with idc check for 1 variable
!     sum_nf23_nbcon         summation for inl,out with nbcon index for 1 variable
!     sum_nf2_m              summation for inl special for udfn_sg_post 2 variables
!
!     grad_press_ls_sum      special summation for grad_press_ls
!     grad_press_ls_c_2d_sum special summation for grad_press_ls_c_2d
!     grad_press_ls__c3d_sum special summation for grad_press_ls_c_3d
!
!------------------------------------------------------------------------------
!
      SUBROUTINE sum_nf(zero,sym, &
                        s1i,s1,   &
                        s2i,s2,   &
                        s3i,s3,   &
                        s4i,s4,   &
                        s5i,s5,   &
                        s6i,s6,   &
                        s7i,s7,   &
                        s8i,s8)
!
      USE Zzone        , ONLY: ncell_fluid
      USE Znum_cell    , ONLY: istart_nb1,                               &
                               ia_nb,icell_nb,right_nb_k,                &
                               nf_number_nb,lens,nf_number_id,istart_nfs
!       
      IMPLICIT NONE
!
!.....Input
      INTEGER,INTENT(IN) :: zero,sym
      REAL(8),INTENT(IN) :: s1i(lens)
      REAL(8),DIMENSION(lens),OPTIONAL :: s2i,s3i,s4i,s5i,s6i,s7i,s8i
!.....Output
      REAL(8),INTENT(OUT) :: s1(ncell_fluid)
      REAL(8),DIMENSION(ncell_fluid),OPTIONAL :: s2,s3,s4,s5,s6,s7,s8
!.....Local variables
      INTEGER :: i,k,nb,nv,ii,nnb,nv0
      INTEGER :: nf_number,len,istart,istart1,i0,i1
      REAL(8) :: sign
      REAL(8) :: s1_s,s2_s,s3_s,s4_s,s5_s,s6_s,s7_s,s8_s
!
      nnb=1
      IF(.not.PRESENT(s2i)) GOTO 100
      nnb=nnb+1
      IF(.not.PRESENT(s3i)) GOTO 100
      nnb=nnb+1
      IF(.not.PRESENT(s4i)) GOTO 100
      nnb=nnb+1
      IF(.not.PRESENT(s5i)) GOTO 100
      nnb=nnb+1
      IF(.not.PRESENT(s6i)) GOTO 100
      nnb=nnb+1
      IF(.not.PRESENT(s7i)) GOTO 100
      nnb=nnb+1
      IF(.not.PRESENT(s8i)) GOTO 100
      nnb=nnb+1
100   CONTINUE
!
      IF(ABS(sym).eq.1) THEN
         nv0=0
      ELSEIF(sym.eq.0) THEN
         nv0=-1
      ENDIF
!
      IF(zero.eq.0) THEN
         IF(nnb.eq.1) THEN
            DO i=1,ncell_fluid
               s1(i)=0.d0
            ENDDO
         ELSEIF(nnb.eq.2) THEN
            DO i=1,ncell_fluid
               s1(i)=0.d0
               s2(i)=0.d0
            ENDDO
         ELSEIF(nnb.eq.3) THEN
            DO i=1,ncell_fluid
               s1(i)=0.d0
               s2(i)=0.d0
               s3(i)=0.d0
            ENDDO
         ELSEIF(nnb.eq.4) THEN
            DO i=1,ncell_fluid
               s1(i)=0.d0
               s2(i)=0.d0
               s3(i)=0.d0
               s4(i)=0.d0
            ENDDO
         ELSEIF(nnb.eq.5) THEN
            DO i=1,ncell_fluid
               s1(i)=0.d0
               s2(i)=0.d0
               s3(i)=0.d0
               s4(i)=0.d0
               s5(i)=0.d0
            ENDDO
         ELSEIF(nnb.eq.6) THEN
            DO i=1,ncell_fluid
               s1(i)=0.d0
               s2(i)=0.d0
               s3(i)=0.d0
               s4(i)=0.d0
               s5(i)=0.d0
               s6(i)=0.d0
            ENDDO
         ELSEIF(nnb.eq.7) THEN
            DO i=1,ncell_fluid
               s1(i)=0.d0
               s2(i)=0.d0
               s3(i)=0.d0
               s4(i)=0.d0
               s5(i)=0.d0
               s6(i)=0.d0
               s7(i)=0.d0
            ENDDO
         ELSEIF(nnb.eq.8) THEN
            DO i=1,ncell_fluid
               s1(i)=0.d0
               s2(i)=0.d0
               s3(i)=0.d0
               s4(i)=0.d0
               s5(i)=0.d0
               s6(i)=0.d0
               s7(i)=0.d0
               s8(i)=0.d0
            ENDDO
         ENDIF
      ENDIF
!
!.....Partial sum for non_k
!
      IF(nf_number_id(0).eq.0) THEN
         IF(abs(sym).eq.1) THEN
            IF(sym.eq.-1) THEN
               sign=-1.d0
            ELSEIF(sym.eq.1) THEN
               sign=1.d0
            ENDIF
         nf_number=-1
         len   =istart_nb1(2,nf_number)
         IF(nnb.eq.1) THEN
            DO nb=1,len
               ii=icell_nb(nb)
               s1_s=0.d0
!DIR$ NOVECTOR
               DO i=ia_nb(nb),ia_nb(nb+1)-1
                  k=right_nb_k(i)
                  s1_s=s1_s+sign*s1i(k)
               ENDDO
               s1(ii)=s1(ii)+s1_s
            ENDDO
         ELSEIF(nnb.eq.2) THEN
            DO nb=1,len
               ii=icell_nb(nb)
               s1_s=0.d0
               s2_s=0.d0
!DIR$ NOVECTOR
               DO i=ia_nb(nb),ia_nb(nb+1)-1
                  k=right_nb_k(i)
                  s1_s=s1_s+sign*s1i(k)
                  s2_s=s2_s+sign*s2i(k)
               ENDDO
               s1(ii)=s1(ii)+s1_s
               s2(ii)=s2(ii)+s2_s
            ENDDO
         ELSEIF(nnb.eq.3) THEN
            DO nb=1,len
               ii=icell_nb(nb)
               s1_s=0.d0
               s2_s=0.d0
               s3_s=0.d0
!DIR$ NOVECTOR
               DO i=ia_nb(nb),ia_nb(nb+1)-1
                  k=right_nb_k(i)
                  s1_s=s1_s+sign*s1i(k)
                  s2_s=s2_s+sign*s2i(k)
                  s3_s=s3_s+sign*s3i(k)
               ENDDO
               s1(ii)=s1(ii)+s1_s
               s2(ii)=s2(ii)+s2_s
               s3(ii)=s3(ii)+s3_s
            ENDDO
         ELSEIF(nnb.eq.4) THEN
            DO nb=1,len
               ii=icell_nb(nb)
               s1_s=0.d0
               s2_s=0.d0
               s3_s=0.d0
               s4_s=0.d0
!DIR$ NOVECTOR
               DO i=ia_nb(nb),ia_nb(nb+1)-1
                  k=right_nb_k(i)
                  s1_s=s1_s+sign*s1i(k)
                  s2_s=s2_s+sign*s2i(k)
                  s3_s=s3_s+sign*s3i(k)
                  s4_s=s4_s+sign*s4i(k)
               ENDDO
               s1(ii)=s1(ii)+s1_s
               s2(ii)=s2(ii)+s2_s
               s3(ii)=s3(ii)+s3_s
               s4(ii)=s4(ii)+s4_s
            ENDDO
         ELSEIF(nnb.eq.5) THEN
            DO nb=1,len
               ii=icell_nb(nb)
               s1_s=0.d0
               s2_s=0.d0
               s3_s=0.d0
               s4_s=0.d0
               s5_s=0.d0
!DIR$ NOVECTOR
               DO i=ia_nb(nb),ia_nb(nb+1)-1
                  k=right_nb_k(i)
                  s1_s=s1_s+sign*s1i(k)
                  s2_s=s2_s+sign*s2i(k)
                  s3_s=s3_s+sign*s3i(k)
                  s4_s=s4_s+sign*s4i(k)
                  s5_s=s5_s+sign*s5i(k)
               ENDDO
               s1(ii)=s1(ii)+s1_s
               s2(ii)=s2(ii)+s2_s
               s3(ii)=s3(ii)+s3_s
               s4(ii)=s4(ii)+s4_s
               s5(ii)=s5(ii)+s5_s
            ENDDO
         ELSEIF(nnb.eq.6) THEN
            DO nb=1,len
               ii=icell_nb(nb)
               s1_s=0.d0
               s2_s=0.d0
               s3_s=0.d0
               s4_s=0.d0
               s5_s=0.d0
               s6_s=0.d0
!DIR$ NOVECTOR
               DO i=ia_nb(nb),ia_nb(nb+1)-1
                  k=right_nb_k(i)
                  s1_s=s1_s+sign*s1i(k)
                  s2_s=s2_s+sign*s2i(k)
                  s3_s=s3_s+sign*s3i(k)
                  s4_s=s4_s+sign*s4i(k)
                  s5_s=s5_s+sign*s5i(k)
                  s6_s=s6_s+sign*s6i(k)
               ENDDO
               s1(ii)=s1(ii)+s1_s
               s2(ii)=s2(ii)+s2_s
               s3(ii)=s3(ii)+s3_s
               s4(ii)=s4(ii)+s4_s
               s5(ii)=s5(ii)+s5_s
               s6(ii)=s6(ii)+s6_s
            ENDDO
         ELSEIF(nnb.eq.7) THEN
            DO nb=1,len
               ii=icell_nb(nb)
               s1_s=0.d0
               s2_s=0.d0
               s3_s=0.d0
               s4_s=0.d0
               s5_s=0.d0
               s6_s=0.d0
               s7_s=0.d0
!DIR$ NOVECTOR
               DO i=ia_nb(nb),ia_nb(nb+1)-1
                  k=right_nb_k(i)
                  s1_s=s1_s+sign*s1i(k)
                  s2_s=s2_s+sign*s2i(k)
                  s3_s=s3_s+sign*s3i(k)
                  s4_s=s4_s+sign*s4i(k)
                  s5_s=s5_s+sign*s5i(k)
                  s6_s=s6_s+sign*s6i(k)
                  s7_s=s7_s+sign*s7i(k)
               ENDDO
               s1(ii)=s1(ii)+s1_s
               s2(ii)=s2(ii)+s2_s
               s3(ii)=s3(ii)+s3_s
               s4(ii)=s4(ii)+s4_s
               s5(ii)=s5(ii)+s5_s
               s6(ii)=s6(ii)+s6_s
               s7(ii)=s7(ii)+s7_s
            ENDDO
         ELSEIF(nnb.eq.8) THEN
            DO nb=1,len
               ii=icell_nb(nb)
               s1_s=0.d0
               s2_s=0.d0
               s3_s=0.d0
               s4_s=0.d0
               s5_s=0.d0
               s6_s=0.d0
               s7_s=0.d0
               s8_s=0.d0
!DIR$ NOVECTOR
               DO i=ia_nb(nb),ia_nb(nb+1)-1
                  k=right_nb_k(i)
                  s1_s=s1_s+sign*s1i(k)
                  s2_s=s2_s+sign*s2i(k)
                  s3_s=s3_s+sign*s3i(k)
                  s4_s=s4_s+sign*s4i(k)
                  s5_s=s5_s+sign*s5i(k)
                  s6_s=s6_s+sign*s6i(k)
                  s7_s=s7_s+sign*s7i(k)
                  s8_s=s8_s+sign*s8i(k)
               ENDDO
               s1(ii)=s1(ii)+s1_s
               s2(ii)=s2(ii)+s2_s
               s3(ii)=s3(ii)+s3_s
               s4(ii)=s4(ii)+s4_s
               s5(ii)=s5(ii)+s5_s
               s6(ii)=s6(ii)+s6_s
               s7(ii)=s7(ii)+s7_s
               s8(ii)=s8(ii)+s8_s
            ENDDO
         ENDIF
      ENDIF
      ENDIF
!
!.....Partial sum for the rest
!
      DO nv=nv0,nf_number_nb
         nf_number=nf_number_id(nv)
         istart =istart_nfs(nv)
         istart1=istart_nb1(1,nf_number)
         len    =istart_nb1(2,nf_number)
         IF(nnb.eq.1) THEN
            DO nb=1,len
               i1=istart1+nb
               ii=icell_nb(i1)
               s1_s=0.d0
!DIR$ NOVECTOR
               DO i=ia_nb(i1),ia_nb(i1+1)-1
                  i0=istart+i
                  s1_s=s1_s+s1i(i0)
               ENDDO
               s1(ii)=s1(ii)+s1_s
            ENDDO
         ELSEIF(nnb.eq.2) THEN
            DO nb=1,len
               i1=istart1+nb
               ii=icell_nb(i1)
               s1_s=0.d0
               s2_s=0.d0
!DIR$ NOVECTOR
               DO i=ia_nb(i1),ia_nb(i1+1)-1
                  i0=istart+i
                  s1_s=s1_s+s1i(i0)
                  s2_s=s2_s+s2i(i0)
               ENDDO
               s1(ii)=s1(ii)+s1_s
               s2(ii)=s2(ii)+s2_s
            ENDDO
         ELSEIF(nnb.eq.3) THEN
            DO nb=1,len
               i1=istart1+nb
               ii=icell_nb(i1)
               s1_s=0.d0
               s2_s=0.d0
               s3_s=0.d0
!DIR$ NOVECTOR
               DO i=ia_nb(i1),ia_nb(i1+1)-1
                  i0=istart+i
                  s1_s=s1_s+s1i(i0)
                  s2_s=s2_s+s2i(i0)
                  s3_s=s3_s+s3i(i0)
               ENDDO
               s1(ii)=s1(ii)+s1_s
               s2(ii)=s2(ii)+s2_s
               s3(ii)=s3(ii)+s3_s
            ENDDO
         ELSEIF(nnb.eq.4) THEN
            DO nb=1,len
               i1=istart1+nb
               ii=icell_nb(i1)
               s1_s=0.d0
               s2_s=0.d0
               s3_s=0.d0
               s4_s=0.d0
!DIR$ NOVECTOR
               DO i=ia_nb(i1),ia_nb(i1+1)-1
                  i0=istart+i
                  s1_s=s1_s+s1i(i0)
                  s2_s=s2_s+s2i(i0)
                  s3_s=s3_s+s3i(i0)
                  s4_s=s4_s+s4i(i0)
               ENDDO
               s1(ii)=s1(ii)+s1_s
               s2(ii)=s2(ii)+s2_s
               s3(ii)=s3(ii)+s3_s
               s4(ii)=s4(ii)+s4_s
            ENDDO
         ELSEIF(nnb.eq.5) THEN
            DO nb=1,len
               i1=istart1+nb
               ii=icell_nb(i1)
               s1_s=0.d0
               s2_s=0.d0
               s3_s=0.d0
               s4_s=0.d0
               s5_s=0.d0
!DIR$ NOVECTOR
               DO i=ia_nb(i1),ia_nb(i1+1)-1
                  i0=istart+i
                  s1_s=s1_s+s1i(i0)
                  s2_s=s2_s+s2i(i0)
                  s3_s=s3_s+s3i(i0)
                  s4_s=s4_s+s4i(i0)
                  s5_s=s5_s+s5i(i0)
               ENDDO
               s1(ii)=s1(ii)+s1_s
               s2(ii)=s2(ii)+s2_s
               s3(ii)=s3(ii)+s3_s
               s4(ii)=s4(ii)+s4_s
               s5(ii)=s5(ii)+s5_s
            ENDDO
         ELSEIF(nnb.eq.6) THEN
            DO nb=1,len
               i1=istart1+nb
               ii=icell_nb(i1)
               s1_s=0.d0
               s2_s=0.d0
               s3_s=0.d0
               s4_s=0.d0
               s5_s=0.d0
               s6_s=0.d0
!DIR$ NOVECTOR
               DO i=ia_nb(i1),ia_nb(i1+1)-1
                  i0=istart+i
                  s1_s=s1_s+s1i(i0)
                  s2_s=s2_s+s2i(i0)
                  s3_s=s3_s+s3i(i0)
                  s4_s=s4_s+s4i(i0)
                  s5_s=s5_s+s5i(i0)
                  s6_s=s6_s+s6i(i0)
               ENDDO
               s1(ii)=s1(ii)+s1_s
               s2(ii)=s2(ii)+s2_s
               s3(ii)=s3(ii)+s3_s
               s4(ii)=s4(ii)+s4_s
               s5(ii)=s5(ii)+s5_s
               s6(ii)=s6(ii)+s6_s
            ENDDO
         ELSEIF(nnb.eq.7) THEN
            DO nb=1,len
               i1=istart1+nb
               ii=icell_nb(i1)
               s1_s=0.d0
               s2_s=0.d0
               s3_s=0.d0
               s4_s=0.d0
               s5_s=0.d0
               s6_s=0.d0
               s7_s=0.d0
!DIR$ NOVECTOR
               DO i=ia_nb(i1),ia_nb(i1+1)-1
                  i0=istart+i
                  s1_s=s1_s+s1i(i0)
                  s2_s=s2_s+s2i(i0)
                  s3_s=s3_s+s3i(i0)
                  s4_s=s4_s+s4i(i0)
                  s5_s=s5_s+s5i(i0)
                  s6_s=s6_s+s6i(i0)
                  s7_s=s7_s+s7i(i0)
               ENDDO
               s1(ii)=s1(ii)+s1_s
               s2(ii)=s2(ii)+s2_s
               s3(ii)=s3(ii)+s3_s
               s4(ii)=s4(ii)+s4_s
               s5(ii)=s5(ii)+s5_s
               s6(ii)=s6(ii)+s6_s
               s7(ii)=s7(ii)+s7_s
            ENDDO
         ELSEIF(nnb.eq.8) THEN
            DO nb=1,len
               i1=istart1+nb
               ii=icell_nb(i1)
               s1_s=0.d0
               s2_s=0.d0
               s3_s=0.d0
               s4_s=0.d0
               s5_s=0.d0
               s6_s=0.d0
               s7_s=0.d0
               s8_s=0.d0
!DIR$ NOVECTOR
               DO i=ia_nb(i1),ia_nb(i1+1)-1
                  i0=istart+i
                  s1_s=s1_s+s1i(i0)
                  s2_s=s2_s+s2i(i0)
                  s3_s=s3_s+s3i(i0)
                  s4_s=s4_s+s4i(i0)
                  s5_s=s5_s+s5i(i0)
                  s6_s=s6_s+s6i(i0)
                  s7_s=s7_s+s7i(i0)
                  s8_s=s8_s+s8i(i0)
               ENDDO
               s1(ii)=s1(ii)+s1_s
               s2(ii)=s2(ii)+s2_s
               s3(ii)=s3(ii)+s3_s
               s4(ii)=s4(ii)+s4_s
               s5(ii)=s5(ii)+s5_s
               s6(ii)=s6(ii)+s6_s
               s7(ii)=s7(ii)+s7_s
               s8(ii)=s8(ii)+s8_s
            ENDDO
         ENDIF
      ENDDO
!
      END SUBROUTINE sum_nf
!
!------------------------------------------------------------------------------
!
      SUBROUTINE sum_nf_ndim(zero,sym,ncell, &
                             s1i,s1,         &
                             s2i,s2,         &
                             s3i,s3,         &
                             s4i,s4)

      USE Zzone        , ONLY: ncell_fluid
      USE Zparam       , ONLY: ndim
      USE Znum_cell    , ONLY: istart_nb1,                                &
                               ia_nb,icell_nb,right_nb_k,                 &
                               nf_number_nb,lens,nf_number_id,istart_nfs
!       
      IMPLICIT NONE
!
!.....Input
      INTEGER,INTENT(IN) :: zero,sym
      INTEGER :: ncell
      REAL(8) :: s1i(lens,ndim)
      REAL(8),DIMENSION(lens,ndim),OPTIONAL :: s2i,s3i,s4i
!.....Output
      REAL(8) :: s1(ncell,ndim)
      REAL(8),DIMENSION(ncell,ndim),OPTIONAL :: s2,s3,s4
!.....Local variables
      INTEGER :: i,k,nb,nv,ii,nnb,nv0
      INTEGER :: nf_number,len,istart,istart1,i0,i1
      REAL(8) :: sign
      REAL(8) :: s1_s1,s1_s2,s1_s3
      REAL(8) :: s2_s1,s2_s2,s2_s3
      REAL(8) :: s3_s1,s3_s2,s3_s3
      REAL(8) :: s4_s1,s4_s2,s4_s3
!
      nnb=1
      IF(.not.PRESENT(s2i)) GOTO 100
      nnb=nnb+1
      IF(.not.PRESENT(s3i)) GOTO 100
      nnb=nnb+1
      IF(.not.PRESENT(s4i)) GOTO 100
      nnb=nnb+1
100   CONTINUE
!
      IF(ABS(sym).eq.1) THEN
         nv0=0
      ELSEIF(sym.eq.0) THEN
         nv0=-1
      ENDIF
!
      IF(ndim.eq.2) THEN
         IF(zero.eq.0) THEN
         IF(nnb.eq.1) THEN
            DO i=1,ncell_fluid
               s1(i,1)=0.d0
               s1(i,2)=0.d0
            ENDDO
         ELSEIF(nnb.eq.2) THEN
            DO i=1,ncell_fluid
               s1(i,1)=0.d0
               s1(i,2)=0.d0
               s2(i,1)=0.d0
               s2(i,2)=0.d0
            ENDDO
         ELSEIF(nnb.eq.3) THEN
            DO i=1,ncell_fluid
               s1(i,1)=0.d0
               s1(i,2)=0.d0
               s2(i,1)=0.d0
               s2(i,2)=0.d0
               s3(i,1)=0.d0
               s3(i,2)=0.d0
            ENDDO
         ELSEIF(nnb.eq.4) THEN
            DO i=1,ncell_fluid
               s1(i,1)=0.d0
               s1(i,2)=0.d0
               s2(i,1)=0.d0
               s2(i,2)=0.d0
               s3(i,1)=0.d0
               s3(i,2)=0.d0
               s4(i,1)=0.d0
               s4(i,2)=0.d0
            ENDDO
         ENDIF
         ENDIF
!
!........Partial sum for non_k
!
         IF(nf_number_id(0).eq.0) THEN
         IF(abs(sym).eq.1) THEN
            IF(sym.eq.-1) THEN
               sign=-1.d0
            ELSEIF(sym.eq.1) THEN
               sign=1.d0
            ENDIF
            nf_number=-1
            len   =istart_nb1(2,nf_number)
            IF(nnb.eq.1) THEN
               DO nb=1,len
                  ii=icell_nb(nb)
                  s1_s1=0.d0
                  s1_s2=0.d0
!DIR$ NOVECTOR
                  DO i=ia_nb(nb),ia_nb(nb+1)-1
                     k=right_nb_k(i)
                     s1_s1=s1_s1+sign*s1i(k,1)
                     s1_s2=s1_s2+sign*s1i(k,2)
                  ENDDO
                  s1(ii,1)=s1(ii,1)+s1_s1
                  s1(ii,2)=s1(ii,2)+s1_s2
               ENDDO
            ELSEIF(nnb.eq.2) THEN
               DO nb=1,len
                  ii=icell_nb(nb)
                  s1_s1=0.d0
                  s1_s2=0.d0
                  s2_s1=0.d0
                  s2_s2=0.d0
                  DO i=ia_nb(nb),ia_nb(nb+1)-1
                     k=right_nb_k(i)
                     s1_s1=s1_s1+sign*s1i(k,1)
                     s1_s2=s1_s2+sign*s1i(k,2)
                     s2_s1=s2_s1+sign*s2i(k,1)
                     s2_s2=s2_s2+sign*s2i(k,2)
                  ENDDO
                  s1(ii,1)=s1(ii,1)+s1_s1
                  s1(ii,2)=s1(ii,2)+s1_s2
                  s2(ii,1)=s2(ii,1)+s2_s1
                  s2(ii,2)=s2(ii,2)+s2_s2
               ENDDO
            ELSEIF(nnb.eq.3) THEN
               DO nb=1,len
                  ii=icell_nb(nb)
                  s1_s1=0.d0
                  s1_s2=0.d0
                  s2_s1=0.d0
                  s2_s2=0.d0
                  s3_s1=0.d0
                  s3_s2=0.d0
!DIR$ NOVECTOR
                  DO i=ia_nb(nb),ia_nb(nb+1)-1
                     k=right_nb_k(i)
                     s1_s1=s1_s1+sign*s1i(k,1)
                     s1_s2=s1_s2+sign*s1i(k,2)
                     s2_s1=s2_s1+sign*s2i(k,1)
                     s2_s2=s2_s2+sign*s2i(k,2)
                     s3_s1=s3_s1+sign*s3i(k,1)
                     s3_s2=s3_s2+sign*s3i(k,2)
                  ENDDO
                  s1(ii,1)=s1(ii,1)+s1_s1
                  s1(ii,2)=s1(ii,2)+s1_s2
                  s2(ii,1)=s2(ii,1)+s2_s1
                  s2(ii,2)=s2(ii,2)+s2_s2
                  s3(ii,1)=s3(ii,1)+s3_s1
                  s3(ii,2)=s3(ii,2)+s3_s2
               ENDDO
            ELSEIF(nnb.eq.4) THEN
               DO nb=1,len
                  ii=icell_nb(nb)
                  s1_s1=0.d0
                  s1_s2=0.d0
                  s2_s1=0.d0
                  s2_s2=0.d0
                  s3_s1=0.d0
                  s3_s2=0.d0
                  s4_s1=0.d0
                  s4_s2=0.d0
!DIR$ NOVECTOR
                  DO i=ia_nb(nb),ia_nb(nb+1)-1
                     k=right_nb_k(i)
                     s1_s1=s1_s1+sign*s1i(k,1)
                     s1_s2=s1_s2+sign*s1i(k,2)
                     s2_s1=s2_s1+sign*s2i(k,1)
                     s2_s2=s2_s2+sign*s2i(k,2)
                     s3_s1=s3_s1+sign*s3i(k,1)
                     s3_s2=s3_s2+sign*s3i(k,2)
                     s4_s1=s4_s1+sign*s4i(k,1)
                     s4_s2=s4_s2+sign*s4i(k,2)
                  ENDDO
                  s1(ii,1)=s1(ii,1)+s1_s1
                  s1(ii,2)=s1(ii,2)+s1_s2
                  s2(ii,1)=s2(ii,1)+s2_s1
                  s2(ii,2)=s2(ii,2)+s2_s2
                  s3(ii,1)=s3(ii,1)+s3_s1
                  s3(ii,2)=s3(ii,2)+s3_s2
                  s4(ii,1)=s4(ii,1)+s4_s1
                  s4(ii,2)=s4(ii,2)+s4_s2
               ENDDO
            ENDIF
         ENDIF
         ENDIF
!
!........Partial sum for the rest
!
         DO nv=nv0,nf_number_nb
            nf_number=nf_number_id(nv)
            istart =istart_nfs(nv)
            istart1=istart_nb1(1,nf_number)
            len    =istart_nb1(2,nf_number)
            IF(nnb.eq.1) THEN
               DO nb=1,len
                  i1=istart1+nb
                  ii=icell_nb(i1)
                  s1_s1=0.d0
                  s1_s2=0.d0
!DIR$ NOVECTOR
                  DO i=ia_nb(i1),ia_nb(i1+1)-1
                     i0=istart+i
                     s1_s1=s1_s1+s1i(i0,1)
                     s1_s2=s1_s2+s1i(i0,2)
                  ENDDO
                  s1(ii,1)=s1(ii,1)+s1_s1
                  s1(ii,2)=s1(ii,2)+s1_s2
               ENDDO
            ELSEIF(nnb.eq.2) THEN
               DO nb=1,len
                  i1=istart1+nb
                  ii=icell_nb(i1)
                  s1_s1=0.d0
                  s1_s2=0.d0
                  s2_s1=0.d0
                  s2_s2=0.d0
!DIR$ NOVECTOR
                  DO i=ia_nb(i1),ia_nb(i1+1)-1
                     i0=istart+i
                     s1_s1=s1_s1+s1i(i0,1)
                     s1_s2=s1_s2+s1i(i0,2)
                     s2_s1=s2_s1+s2i(i0,1)
                     s2_s2=s2_s2+s2i(i0,2)
                  ENDDO
                  s1(ii,1)=s1(ii,1)+s1_s1
                  s1(ii,2)=s1(ii,2)+s1_s2
                  s2(ii,1)=s2(ii,1)+s2_s1
                  s2(ii,2)=s2(ii,2)+s2_s2
               ENDDO
            ELSEIF(nnb.eq.3) THEN
               DO nb=1,len
                  i1=istart1+nb
                  ii=icell_nb(i1)
                  s1_s1=0.d0
                  s1_s2=0.d0
                  s2_s1=0.d0
                  s2_s2=0.d0
                  s3_s1=0.d0
                  s3_s2=0.d0
!DIR$ NOVECTOR
                  DO i=ia_nb(i1),ia_nb(i1+1)-1
                     i0=istart+i
                     s1_s1=s1_s1+s1i(i0,1)
                     s1_s2=s1_s2+s1i(i0,2)
                     s2_s1=s2_s1+s2i(i0,1)
                     s2_s2=s2_s2+s2i(i0,2)
                     s3_s1=s3_s1+s3i(i0,1)
                     s3_s2=s3_s2+s3i(i0,2)
                  ENDDO
                  s1(ii,1)=s1(ii,1)+s1_s1
                  s1(ii,2)=s1(ii,2)+s1_s2
                  s2(ii,1)=s2(ii,1)+s2_s1
                  s2(ii,2)=s2(ii,2)+s2_s2
                  s3(ii,1)=s3(ii,1)+s3_s1
                  s3(ii,2)=s3(ii,2)+s3_s2
               ENDDO
            ELSEIF(nnb.eq.4) THEN
               DO nb=1,len
                  i1=istart1+nb
                  ii=icell_nb(i1)
                  s1_s1=0.d0
                  s1_s2=0.d0
                  s2_s1=0.d0
                  s2_s2=0.d0
                  s3_s1=0.d0
                  s3_s2=0.d0
                  s4_s1=0.d0
                  s4_s2=0.d0
!DIR$ NOVECTOR
                  DO i=ia_nb(i1),ia_nb(i1+1)-1
                     i0=istart+i
                     s1_s1=s1_s1+s1i(i0,1)
                     s1_s2=s1_s2+s1i(i0,2)
                     s2_s1=s2_s1+s2i(i0,1)
                     s2_s2=s2_s2+s2i(i0,2)
                     s3_s1=s3_s1+s3i(i0,1)
                     s3_s2=s3_s2+s3i(i0,2)
                     s4_s1=s4_s1+s4i(i0,1)
                     s4_s2=s4_s2+s4i(i0,2)
                  ENDDO
                  s1(ii,1)=s1(ii,1)+s1_s1
                  s1(ii,2)=s1(ii,2)+s1_s2
                  s2(ii,1)=s2(ii,1)+s2_s1
                  s2(ii,2)=s2(ii,2)+s2_s2
                  s3(ii,1)=s3(ii,1)+s3_s1
                  s3(ii,2)=s3(ii,2)+s3_s2
                  s4(ii,1)=s4(ii,1)+s4_s1
                  s4(ii,2)=s4(ii,2)+s4_s2
               ENDDO
            ENDIF
         ENDDO
      ELSE
         IF(zero.eq.0) THEN
         IF(nnb.eq.1) THEN
            DO i=1,ncell_fluid
               s1(i,1)=0.d0
               s1(i,2)=0.d0
               s1(i,3)=0.d0
            ENDDO
         ELSEIF(nnb.eq.2) THEN
            DO i=1,ncell_fluid
               s1(i,1)=0.d0
               s1(i,2)=0.d0
               s1(i,3)=0.d0
               s2(i,1)=0.d0
               s2(i,2)=0.d0
               s2(i,3)=0.d0
            ENDDO
         ELSEIF(nnb.eq.3) THEN
            DO i=1,ncell_fluid
               s1(i,1)=0.d0
               s1(i,2)=0.d0
               s1(i,3)=0.d0
               s2(i,1)=0.d0
               s2(i,2)=0.d0
               s2(i,3)=0.d0
               s3(i,1)=0.d0
               s3(i,2)=0.d0
               s3(i,3)=0.d0
            ENDDO
         ELSEIF(nnb.eq.4) THEN
            DO i=1,ncell_fluid
               s1(i,1)=0.d0
               s1(i,2)=0.d0
               s1(i,3)=0.d0
               s2(i,1)=0.d0
               s2(i,2)=0.d0
               s2(i,3)=0.d0
               s3(i,1)=0.d0
               s3(i,2)=0.d0
               s3(i,3)=0.d0
               s4(i,1)=0.d0
               s4(i,2)=0.d0
               s4(i,3)=0.d0
            ENDDO
         ENDIF
         ENDIF
!
!........Partial sum for non_k
!
         IF(nf_number_id(0).eq.0) THEN
         IF(abs(sym).eq.1) THEN
            IF(sym.eq.-1) THEN
               sign=-1.d0
            ELSEIF(sym.eq.1) THEN
               sign=1.d0
            ENDIF
            nf_number=-1
            len   =istart_nb1(2,nf_number)
            IF(nnb.eq.1) THEN
               DO nb=1,len
                  ii=icell_nb(nb)
                  s1_s1=0.d0
                  s1_s2=0.d0
                  s1_s3=0.d0
!DIR$ NOVECTOR
                  DO i=ia_nb(nb),ia_nb(nb+1)-1
                     k=right_nb_k(i)
                     s1_s1=s1_s1+sign*s1i(k,1)
                     s1_s2=s1_s2+sign*s1i(k,2)
                     s1_s3=s1_s3+sign*s1i(k,3)
                  ENDDO
                  s1(ii,1)=s1(ii,1)+s1_s1
                  s1(ii,2)=s1(ii,2)+s1_s2
                  s1(ii,3)=s1(ii,3)+s1_s3
               ENDDO
            ELSEIF(nnb.eq.2) THEN
               DO nb=1,len
                  ii=icell_nb(nb)
                  s1_s1=0.d0
                  s1_s2=0.d0
                  s1_s3=0.d0
                  s2_s1=0.d0
                  s2_s2=0.d0
                  s2_s3=0.d0
!DIR$ NOVECTOR
                  DO i=ia_nb(nb),ia_nb(nb+1)-1
                     k=right_nb_k(i)
                     s1_s1=s1_s1+sign*s1i(k,1)
                     s1_s2=s1_s2+sign*s1i(k,2)
                     s1_s3=s1_s3+sign*s1i(k,3)
                     s2_s1=s2_s1+sign*s2i(k,1)
                     s2_s2=s2_s2+sign*s2i(k,2)
                     s2_s3=s2_s3+sign*s2i(k,3)
                  ENDDO
                  s1(ii,1)=s1(ii,1)+s1_s1
                  s1(ii,2)=s1(ii,2)+s1_s2
                  s1(ii,3)=s1(ii,3)+s1_s3
                  s2(ii,1)=s2(ii,1)+s2_s1
                  s2(ii,2)=s2(ii,2)+s2_s2
                  s2(ii,3)=s2(ii,3)+s2_s3
               ENDDO
            ELSEIF(nnb.eq.3) THEN
               DO nb=1,len
                  ii=icell_nb(nb)
                  s1_s1=0.d0
                  s1_s2=0.d0
                  s1_s3=0.d0
                  s2_s1=0.d0
                  s2_s2=0.d0
                  s2_s3=0.d0
                  s3_s1=0.d0
                  s3_s2=0.d0
                  s3_s3=0.d0
!DIR$ NOVECTOR
                  DO i=ia_nb(nb),ia_nb(nb+1)-1
                     k=right_nb_k(i)
                     s1_s1=s1_s1+sign*s1i(k,1)
                     s1_s2=s1_s2+sign*s1i(k,2)
                     s1_s3=s1_s3+sign*s1i(k,3)
                     s2_s1=s2_s1+sign*s2i(k,1)
                     s2_s2=s2_s2+sign*s2i(k,2)
                     s2_s3=s2_s3+sign*s2i(k,3)
                     s3_s1=s3_s1+sign*s3i(k,1)
                     s3_s2=s3_s2+sign*s3i(k,2)
                     s3_s3=s3_s3+sign*s3i(k,3)
                  ENDDO
                  s1(ii,1)=s1(ii,1)+s1_s1
                  s1(ii,2)=s1(ii,2)+s1_s2
                  s1(ii,3)=s1(ii,3)+s1_s3
                  s2(ii,1)=s2(ii,1)+s2_s1
                  s2(ii,2)=s2(ii,2)+s2_s2
                  s2(ii,3)=s2(ii,3)+s2_s3
                  s3(ii,1)=s3(ii,1)+s3_s1
                  s3(ii,2)=s3(ii,2)+s3_s2
                  s3(ii,3)=s3(ii,3)+s3_s3
               ENDDO
            ELSEIF(nnb.eq.4) THEN
               DO nb=1,len
                  ii=icell_nb(nb)
                  s1_s1=0.d0
                  s1_s2=0.d0
                  s1_s3=0.d0
                  s2_s1=0.d0
                  s2_s2=0.d0
                  s2_s3=0.d0
                  s3_s1=0.d0
                  s3_s2=0.d0
                  s3_s3=0.d0
                  s4_s1=0.d0
                  s4_s2=0.d0
                  s4_s3=0.d0
!DIR$ NOVECTOR
                  DO i=ia_nb(nb),ia_nb(nb+1)-1
                     k=right_nb_k(i)
                     s1_s1=s1_s1+sign*s1i(k,1)
                     s1_s2=s1_s2+sign*s1i(k,2)
                     s1_s3=s1_s3+sign*s1i(k,3)
                     s2_s1=s2_s1+sign*s2i(k,1)
                     s2_s2=s2_s2+sign*s2i(k,2)
                     s2_s3=s2_s3+sign*s2i(k,3)
                     s3_s1=s3_s1+sign*s3i(k,1)
                     s3_s2=s3_s2+sign*s3i(k,2)
                     s3_s3=s3_s3+sign*s3i(k,3)
                     s4_s1=s4_s1+sign*s4i(k,1)
                     s4_s2=s4_s2+sign*s4i(k,2)
                     s4_s3=s4_s3+sign*s4i(k,3)
                  ENDDO
                  s1(ii,1)=s1(ii,1)+s1_s1
                  s1(ii,2)=s1(ii,2)+s1_s2
                  s1(ii,3)=s1(ii,3)+s1_s3
                  s2(ii,1)=s2(ii,1)+s2_s1
                  s2(ii,2)=s2(ii,2)+s2_s2
                  s2(ii,3)=s2(ii,3)+s2_s3
                  s3(ii,1)=s3(ii,1)+s3_s1
                  s3(ii,2)=s3(ii,2)+s3_s2
                  s3(ii,3)=s3(ii,3)+s3_s3
                  s4(ii,1)=s4(ii,1)+s4_s1
                  s4(ii,2)=s4(ii,2)+s4_s2
                  s4(ii,3)=s4(ii,3)+s4_s3
               ENDDO
            ENDIF
         ENDIF
         ENDIF
!
!........Partial sum for the rest
!
         DO nv=nv0,nf_number_nb
            nf_number=nf_number_id(nv)
            istart =istart_nfs(nv)
            istart1=istart_nb1(1,nf_number)
            len    =istart_nb1(2,nf_number)
            IF(nnb.eq.1) THEN
               DO nb=1,len
                  i1=istart1+nb
                  ii=icell_nb(i1)
                  s1_s1=0.d0
                  s1_s2=0.d0
                  s1_s3=0.d0
!DIR$ NOVECTOR
                  DO i=ia_nb(i1),ia_nb(i1+1)-1
                     i0=istart+i
                     s1_s1=s1_s1+s1i(i0,1)
                     s1_s2=s1_s2+s1i(i0,2)
                     s1_s3=s1_s3+s1i(i0,3)
                  ENDDO
                  s1(ii,1)=s1(ii,1)+s1_s1
                  s1(ii,2)=s1(ii,2)+s1_s2
                  s1(ii,3)=s1(ii,3)+s1_s3
               ENDDO
            ELSEIF(nnb.eq.2) THEN
               DO nb=1,len
                  i1=istart1+nb
                  ii=icell_nb(i1)
                  s1_s1=0.d0
                  s1_s2=0.d0
                  s1_s3=0.d0
                  s2_s1=0.d0
                  s2_s2=0.d0
                  s2_s3=0.d0
!DIR$ NOVECTOR
                  DO i=ia_nb(i1),ia_nb(i1+1)-1
                     i0=istart+i
                     s1_s1=s1_s1+s1i(i0,1)
                     s1_s2=s1_s2+s1i(i0,2)
                     s1_s3=s1_s3+s1i(i0,3)
                     s2_s1=s2_s1+s2i(i0,1)
                     s2_s2=s2_s2+s2i(i0,2)
                     s2_s3=s2_s3+s2i(i0,3)
                  ENDDO
                  s1(ii,1)=s1(ii,1)+s1_s1
                  s1(ii,2)=s1(ii,2)+s1_s2
                  s1(ii,3)=s1(ii,3)+s1_s3
                  s2(ii,1)=s2(ii,1)+s2_s1
                  s2(ii,2)=s2(ii,2)+s2_s2
                  s2(ii,3)=s2(ii,3)+s2_s3
               ENDDO
            ELSEIF(nnb.eq.3) THEN
               DO nb=1,len
                  i1=istart1+nb
                  ii=icell_nb(i1)
                  s1_s1=0.d0
                  s1_s2=0.d0
                  s1_s3=0.d0
                  s2_s1=0.d0
                  s2_s2=0.d0
                  s2_s3=0.d0
                  s3_s1=0.d0
                  s3_s2=0.d0
                  s3_s3=0.d0
!DIR$ NOVECTOR
                  DO i=ia_nb(i1),ia_nb(i1+1)-1
                     i0=istart+i
                     s1_s1=s1_s1+s1i(i0,1)
                     s1_s2=s1_s2+s1i(i0,2)
                     s1_s3=s1_s3+s1i(i0,3)
                     s2_s1=s2_s1+s2i(i0,1)
                     s2_s2=s2_s2+s2i(i0,2)
                     s2_s3=s2_s3+s2i(i0,3)
                     s3_s1=s3_s1+s3i(i0,1)
                     s3_s2=s3_s2+s3i(i0,2)
                     s3_s3=s3_s3+s3i(i0,3)
                  ENDDO
                  s1(ii,1)=s1(ii,1)+s1_s1
                  s1(ii,2)=s1(ii,2)+s1_s2
                  s1(ii,3)=s1(ii,3)+s1_s3
                  s2(ii,1)=s2(ii,1)+s2_s1
                  s2(ii,2)=s2(ii,2)+s2_s2
                  s2(ii,3)=s2(ii,3)+s2_s3
                  s3(ii,1)=s3(ii,1)+s3_s1
                  s3(ii,2)=s3(ii,2)+s3_s2
                  s3(ii,3)=s3(ii,3)+s3_s3
               ENDDO
            ELSEIF(nnb.eq.4) THEN
               DO nb=1,len
                  i1=istart1+nb
                  ii=icell_nb(i1)
                  s1_s1=0.d0
                  s1_s2=0.d0
                  s1_s3=0.d0
                  s2_s1=0.d0
                  s2_s2=0.d0
                  s2_s3=0.d0
                  s3_s1=0.d0
                  s3_s2=0.d0
                  s3_s3=0.d0
                  s4_s1=0.d0
                  s4_s2=0.d0
                  s4_s3=0.d0
!DIR$ NOVECTOR
                  DO i=ia_nb(i1),ia_nb(i1+1)-1
                     i0=istart+i
                     s1_s1=s1_s1+s1i(i0,1)
                     s1_s2=s1_s2+s1i(i0,2)
                     s1_s3=s1_s3+s1i(i0,3)
                     s2_s1=s2_s1+s2i(i0,1)
                     s2_s2=s2_s2+s2i(i0,2)
                     s2_s3=s2_s3+s2i(i0,3)
                     s3_s1=s3_s1+s3i(i0,1)
                     s3_s2=s3_s2+s3i(i0,2)
                     s3_s3=s3_s3+s3i(i0,3)
                     s4_s1=s4_s1+s4i(i0,1)
                     s4_s2=s4_s2+s4i(i0,2)
                     s4_s3=s4_s3+s4i(i0,3)
                  ENDDO
                  s1(ii,1)=s1(ii,1)+s1_s1
                  s1(ii,2)=s1(ii,2)+s1_s2
                  s1(ii,3)=s1(ii,3)+s1_s3
                  s2(ii,1)=s2(ii,1)+s2_s1
                  s2(ii,2)=s2(ii,2)+s2_s2
                  s2(ii,3)=s2(ii,3)+s2_s3
                  s3(ii,1)=s3(ii,1)+s3_s1
                  s3(ii,2)=s3(ii,2)+s3_s2
                  s3(ii,3)=s3(ii,3)+s3_s3
                  s4(ii,1)=s4(ii,1)+s4_s1
                  s4(ii,2)=s4(ii,2)+s4_s2
                  s4(ii,3)=s4(ii,3)+s4_s3
               ENDDO
            ENDIF
         ENDDO
      ENDIF
!
      END SUBROUTINE sum_nf_ndim
!
!------------------------------------------------------------------------------
!
      SUBROUTINE sum_nf_ndim2(zero,sym,ncell, &
                              s1i,s1)
!
      USE Zzone        , ONLY: ncell_fluid
      USE Zparam       , ONLY: ndim
      USE Znum_cell    , ONLY: istart_nb1,                                &
                               ia_nb,icell_nb,right_nb_k,                 &
                               nf_number_nb,lens,nf_number_id,istart_nfs
!       
      IMPLICIT NONE
!
!.....Input
      INTEGER,INTENT(IN) :: zero,sym
      INTEGER :: ncell
      REAL(8) :: s1i(lens,ndim,ndim)
!.....Output
      REAL(8) :: s1(ncell,ndim,ndim)
!.....Local variables
      INTEGER :: i,k,nb,nv,ii,nv0
      INTEGER :: nf_number,len,istart,istart1,i0,i1
      REAL(8) :: sign
      REAL(8) :: s1_s11,s1_s21,s1_s31
      REAL(8) :: s1_s12,s1_s22,s1_s32
      REAL(8) :: s1_s13,s1_s23,s1_s33
!
      IF(abs(sym).eq.1) THEN
         nv0=0
      ELSEIF(sym.eq.0) THEN
         nv0=-1
      ENDIF
!
      IF(ndim.eq.2) THEN
         IF(zero.eq.0) THEN
         DO i=1,ncell_fluid
            s1(i,1,1)=0.d0
            s1(i,2,1)=0.d0
            s1(i,1,2)=0.d0
            s1(i,2,2)=0.d0
         ENDDO
         ENDIF
!
!........Partial sum for non_k
!
         IF(nf_number_id(0).eq.0) THEN
         IF(abs(sym).eq.1) THEN
            IF(sym.eq.-1) THEN
               sign=-1.d0
            ELSEIF(sym.eq.1) THEN
               sign=1.d0
            ENDIF
            nf_number=-1
            len   =istart_nb1(2,nf_number)
            DO nb=1,len
               ii=icell_nb(nb)
               s1_s11=0.d0
               s1_s21=0.d0
               s1_s12=0.d0
               s1_s22=0.d0
!DIR$ NOVECTOR
               DO i=ia_nb(nb),ia_nb(nb+1)-1
                  k=right_nb_k(i)
                  s1_s11=s1_s11+sign*s1i(k,1,1)
                  s1_s21=s1_s21+sign*s1i(k,2,1)
                  s1_s12=s1_s12+sign*s1i(k,1,2)
                  s1_s22=s1_s22+sign*s1i(k,2,2)
               ENDDO
               s1(ii,1,1)=s1(ii,1,1)+s1_s11
               s1(ii,2,1)=s1(ii,2,1)+s1_s21
               s1(ii,1,2)=s1(ii,1,2)+s1_s12
               s1(ii,2,2)=s1(ii,2,2)+s1_s22
            ENDDO
         ENDIF
         ENDIF
!
!........Partial sum for the rest
!
         DO nv=nv0,nf_number_nb
            nf_number=nf_number_id(nv)
            istart =istart_nfs(nv)
            istart1=istart_nb1(1,nf_number)
            len    =istart_nb1(2,nf_number)
            DO nb=1,len
               i1=istart1+nb
               ii=icell_nb(i1)
               s1_s11=0.d0
               s1_s21=0.d0
               s1_s12=0.d0
               s1_s22=0.d0
!DIR$ NOVECTOR
               DO i=ia_nb(i1),ia_nb(i1+1)-1
                  i0=istart+i
                  s1_s11=s1_s11+s1i(i0,1,1)
                  s1_s21=s1_s21+s1i(i0,2,1)
                  s1_s12=s1_s12+s1i(i0,1,2)
                  s1_s22=s1_s22+s1i(i0,2,2)
               ENDDO
               s1(ii,1,1)=s1(ii,1,1)+s1_s11
               s1(ii,2,1)=s1(ii,2,1)+s1_s21
               s1(ii,1,2)=s1(ii,1,2)+s1_s12
               s1(ii,2,2)=s1(ii,2,2)+s1_s22
            ENDDO
         ENDDO
      ELSE
         IF(zero.eq.0) THEN
         DO i=1,ncell_fluid
            s1(i,1,1)=0.d0
            s1(i,2,1)=0.d0
            s1(i,3,1)=0.d0
            s1(i,1,2)=0.d0
            s1(i,2,2)=0.d0
            s1(i,3,2)=0.d0
            s1(i,1,3)=0.d0
            s1(i,2,3)=0.d0
            s1(i,3,3)=0.d0
         ENDDO
         ENDIF
!
!........Partial sum for non_k
!
         IF(nf_number_id(0).eq.0) THEN
         IF(abs(sym).eq.1) THEN
            IF(sym.eq.-1) THEN
               sign=-1.d0
            ELSEIF(sym.eq.1) THEN
               sign=1.d0
            ENDIF
            nf_number=-1
            len   =istart_nb1(2,nf_number)
            DO nb=1,len
               ii=icell_nb(nb)
               s1_s11=0.d0
               s1_s21=0.d0
               s1_s31=0.d0
               s1_s12=0.d0
               s1_s22=0.d0
               s1_s32=0.d0
               s1_s13=0.d0
               s1_s23=0.d0
               s1_s33=0.d0
!DIR$ NOVECTOR
               DO i=ia_nb(nb),ia_nb(nb+1)-1
                  k=right_nb_k(i)
                  s1_s11=s1_s11+sign*s1i(k,1,1)
                  s1_s21=s1_s21+sign*s1i(k,2,1)
                  s1_s31=s1_s31+sign*s1i(k,3,1)
                  s1_s12=s1_s12+sign*s1i(k,1,2)
                  s1_s22=s1_s22+sign*s1i(k,2,2)
                  s1_s32=s1_s32+sign*s1i(k,3,2)
                  s1_s13=s1_s13+sign*s1i(k,1,3)
                  s1_s23=s1_s23+sign*s1i(k,2,3)
                  s1_s33=s1_s33+sign*s1i(k,3,3)
               ENDDO
               s1(ii,1,1)=s1(ii,1,1)+s1_s11
               s1(ii,2,1)=s1(ii,2,1)+s1_s21
               s1(ii,3,1)=s1(ii,3,1)+s1_s31
               s1(ii,1,2)=s1(ii,1,2)+s1_s12
               s1(ii,2,2)=s1(ii,2,2)+s1_s22
               s1(ii,3,2)=s1(ii,3,2)+s1_s32
               s1(ii,1,3)=s1(ii,1,3)+s1_s13
               s1(ii,2,3)=s1(ii,2,3)+s1_s23
               s1(ii,3,3)=s1(ii,3,3)+s1_s33
            ENDDO
         ENDIF
         ENDIF
!
!........Partial sum for the rest
!
         DO nv=nv0,nf_number_nb
            nf_number=nf_number_id(nv)
            istart =istart_nfs(nv)
            istart1=istart_nb1(1,nf_number)
            len    =istart_nb1(2,nf_number)
            DO nb=1,len
               i1=istart1+nb
               ii=icell_nb(i1)
               s1_s11=0.d0
               s1_s21=0.d0
               s1_s31=0.d0
               s1_s12=0.d0
               s1_s22=0.d0
               s1_s32=0.d0
               s1_s13=0.d0
               s1_s23=0.d0
               s1_s33=0.d0
!DIR$ NOVECTOR
               DO i=ia_nb(i1),ia_nb(i1+1)-1
                  i0=istart+i
                  s1_s11=s1_s11+s1i(i0,1,1)
                  s1_s21=s1_s21+s1i(i0,2,1)
                  s1_s31=s1_s31+s1i(i0,3,1)
                  s1_s12=s1_s12+s1i(i0,1,2)
                  s1_s22=s1_s22+s1i(i0,2,2)
                  s1_s32=s1_s32+s1i(i0,3,2)
                  s1_s13=s1_s13+s1i(i0,1,3)
                  s1_s23=s1_s23+s1i(i0,2,3)
                  s1_s33=s1_s33+s1i(i0,3,3)
               ENDDO
               s1(ii,1,1)=s1(ii,1,1)+s1_s11
               s1(ii,2,1)=s1(ii,2,1)+s1_s21
               s1(ii,3,1)=s1(ii,3,1)+s1_s31
               s1(ii,1,2)=s1(ii,1,2)+s1_s12
               s1(ii,2,2)=s1(ii,2,2)+s1_s22
               s1(ii,3,2)=s1(ii,3,2)+s1_s32
               s1(ii,1,3)=s1(ii,1,3)+s1_s13
               s1(ii,2,3)=s1(ii,2,3)+s1_s23
               s1(ii,3,3)=s1(ii,3,3)+s1_s33
            ENDDO
         ENDDO
      ENDIF
!
      END SUBROUTINE sum_nf_ndim2
!
!------------------------------------------------------------------------------
!
      SUBROUTINE sum_nf_solid(isign,zero,sym, &
                              s1i,s1)
!
      USE Zzone        , ONLY: ncell_cond
      USE Znum_cell    , ONLY: istartc_nb1,                              &
                               iac_nb,icellc_nb,                         &
                               nf_number_nb,lens,nf_number_id,istart_nfs
!
      IMPLICIT NONE
!
!.....Input
      INTEGER,INTENT(IN) :: isign,zero,sym
      REAL(8),DIMENSION(lens) :: s1i
!.....Output
      REAL(8),DIMENSION(ncell_cond) :: s1
!.....Local variables
      INTEGER :: i,ii,nb,nv0
      INTEGER :: nv,nf_number,len,istart,istart1,i0,i1
      REAL(8) :: sign
      REAL(8) :: s1_s
!
      IF(sym.eq.1) THEN
         nv0=1
      ELSEIF(sym.eq.0) THEN
         nv0=-1
      ENDIF
!
      IF(zero.eq.0) THEN
      DO i=1,ncell_cond
         s1(i)=0.d0
      ENDDO
      ENDIF
      IF(isign.eq.-1) THEN
         sign=-1.d0
      ELSEIF(isign.eq.1) THEN
         sign=1.d0
      ENDIF
!
!.....Partial sum for the rest
!
      DO nv=nv0,nf_number_nb
         nf_number=nf_number_id(nv)
         istart =istart_nfs(nv)
         istart1=istartc_nb1(1,nf_number)
         len    =istartc_nb1(2,nf_number)
         DO nb=1,len
            i1=istart1+nb
            ii=icellc_nb(i1)
            s1_s=0.d0
!DIR$ NOVECTOR
            DO i=iac_nb(i1),iac_nb(i1+1)-1
               i0=istart+i
               s1_s=s1_s+s1i(i0)
            ENDDO
            s1(ii)=s1(ii)+s1_s
         ENDDO
      ENDDO
!
      END SUBROUTINE sum_nf_solid
!
!------------------------------------------------------------------------------
!
      SUBROUTINE max_nf(s1i,s1)
!
      USE Zzone         , ONLY: ncell_fluid,icore
      USE Zvec_param    , ONLY: nf_tot
      USE Znum_cell     , ONLY: istart_nf,istart_nb1,     &
                                ia_nb,icell_nb,right_nb_k
      USE Zcoord3       , ONLY: porosity
!
      IMPLICIT NONE
!.....Input
      REAL(8) :: s1i(nf_tot)
!.....Output
      REAL(8) :: s1(ncell_fluid)
!.....Local variables
      INTEGER :: i,k,nb
      INTEGER :: ii
      INTEGER :: nf_number,istart,len,istart1,i0,i1
      REAL(8) :: s1_s
!
      DO i=1,ncell_fluid
         s1(i)=0.d0
      ENDDO
!
!.....Partial max for non_k
!
      nf_number=-1
      len   =istart_nb1(2,nf_number)
      DO nb=1,len
         ii=icell_nb(nb)
         IF(icore(ii).eq.1 .or. porosity(ii).ge.1.0d0) cycle
         s1_s=0.d0
!DIR$ NOVECTOR
         DO i=ia_nb(nb),ia_nb(nb+1)-1
            k=right_nb_k(i)
            s1_s=max(s1_s,s1i(k))
         ENDDO
         s1(ii)=s1_s
      ENDDO
!
!.....Partial max for the rest
!
      DO nf_number=0,8
         istart =istart_nf(1,nf_number)
         istart1=istart_nb1(1,nf_number)
         len    =istart_nb1(2,nf_number)
         DO nb=1,len
            i1=istart1+nb
            ii=icell_nb(i1)
            IF(icore(ii).eq.1 .or. porosity(ii).ge.1.0d0) cycle
            s1_s=0.d0
!DIR$ NOVECTOR
            DO i=ia_nb(i1),ia_nb(i1+1)-1
               i0=istart+i
               s1_s=max(s1_s,s1i(i0))
            ENDDO
            s1(ii)=max(s1(ii),s1_s)
         ENDDO
      ENDDO
!
      END SUBROUTINE max_nf
!
!------------------------------------------------------------------------------
!
      SUBROUTINE sum_nf0_idc(s1i,s1,idc_c, &
                            s2i,s2,idc_h)
!
      USE Zmpi          , ONLY: ncell_fp
      USE Zzone         , ONLY: ncell_fluid
      USE Znum_cell     , ONLY: right_nb_k,istart_nb1,       &
                                ia_nb,icell_nb,              &
                                lens,nf_number_id,istart_nfs
      USE Zvec_index    , ONLY: left_nf,right_non
      USE Zsg           , ONLY: idc
!
      IMPLICIT NONE
!.....Input
      REAL(8),DIMENSION(lens) :: s1i,s2i
      INTEGER,DIMENSION(ncell_fp) :: idc_c,idc_h
!.....Output
      REAL(8),DIMENSION(ncell_fluid) :: s1,s2
!.....Local variables
      INTEGER :: i,k,nb
      INTEGER :: ii,kk
      INTEGER :: nv,nf_number,len,istart,istart1,i0,i1
      REAL(8) :: s1_s,s2_s
!
      DO i=1,ncell_fluid
         s1(i)=0.d0
         s2(i)=0.d0
      ENDDO
!
!.....Partial sum diagonal and source for non_k
!
      nv=-1
      nf_number=nf_number_id(nv)
      len   =istart_nb1(2,nf_number)
      DO nb=1,len
         ii=icell_nb(nb)
         IF(idc_c(ii).gt.0) THEN
         s1_s=0.d0
!DIR$ NOVECTOR
         DO i=ia_nb(nb),ia_nb(nb+1)-1
            k=right_nb_k(i)
            kk=left_nf(k)
            IF(idc(kk).eq.1) s1_s=s1_s-s1i(i)
         ENDDO
         s1(ii)=s1_s
         ENDIF
         IF(idc_h(ii).gt.0) THEN
         s2_s=0.d0
!DIR$ NOVECTOR
         DO i=ia_nb(nb),ia_nb(nb+1)-1
            k=right_nb_k(i)
            kk=left_nf(k)
            IF(idc(kk).eq.1) s2_s=s2_s-s2i(i)
         ENDDO
         s2(ii)=s2_s
         ENDIF
      ENDDO
!
!.....Partial sum diagonal only for non_i
!
      nv=0
      nf_number=nf_number_id(nv)
      istart =istart_nfs(nv)
      istart1=istart_nb1(1,nf_number)
      len    =istart_nb1(2,nf_number)
      DO nb=1,len
         i1=istart1+nb
         ii=icell_nb(i1)
         IF(idc_c(ii).gt.0) THEN
         s1_s=0.d0
!DIR$ NOVECTOR
         DO i=ia_nb(i1),ia_nb(i1+1)-1
            i0=istart+i
            kk=right_non(i)
            IF(idc(kk).eq.1) s1_s=s1_s+s1i(i0)
         ENDDO
         s1(ii)=s1(ii)+s1_s
         ENDIF
         IF(idc_h(ii).gt.0) THEN
         s2_s=0.d0
!DIR$ NOVECTOR
         DO i=ia_nb(i1),ia_nb(i1+1)-1
            i0=istart+i
            kk=right_non(i)
            IF(idc(kk).eq.1) s2_s=s2_s+s2i(i0)
         ENDDO
         s2(ii)=s2(ii)+s2_s
         ENDIF
      ENDDO
!
      END SUBROUTINE sum_nf0_idc
!
!------------------------------------------------------------------------------
!
      SUBROUTINE sum_nf2_m(s1i,s1, &
                           s2i,s2)
!
      USE Zzone         , ONLY: ncell_fluid
      USE Znum_cell     , ONLY: istart_nb1,                               &
                                ia_nb,icell_nb,                           &
                                istart_nbcon_nf,                          &
                                nf_number_nb,lens,nf_number_id,istart_nfs
      USE Zvec_index    , ONLY: nbcon_nf
!
      IMPLICIT NONE
!.....Input
      REAL(8),DIMENSION(lens) :: s1i,s2i
!.....Output
      REAL(8) :: s1(ncell_fluid),s2(ncell_fluid)
!.....Local variables
      INTEGER :: i,k,nb
      INTEGER :: ii
      INTEGER :: nv,nf_number,len,istart,istart1,i0,i1
      INTEGER :: istart3,i3
      REAL(8) :: s1_s,s2_s
!
      DO i=1,ncell_fluid
         s1(i)=0.d0
         s2(i)=0.d0
      ENDDO
!
!.....Partial sum for all
!
      DO nv=0,nf_number_nb
         nf_number=nf_number_id(nv)
         istart =istart_nfs(nv)
         istart1=istart_nb1(1,nf_number)
         len    =istart_nb1(2,nf_number)
         istart3=istart_nbcon_nf(nf_number)
         DO nb=1,len
            i1=istart1+nb
            ii=icell_nb(i1)
            s1_s=0.d0
            s2_s=0.d0
!DIR$ NOVECTOR
            DO i=ia_nb(i1),ia_nb(i1+1)-1
               i0=istart+i
               i3=istart3+i
               k=nbcon_nf(i3)
               IF(k.eq.1)THEN
                  s1_s=s1_s-s1i(i0)
               ELSEIF (k.eq.2)THEN
                  s2_s=s2_s-s2i(i0)
               ENDIF
            ENDDO
            s1(ii)=s1(ii)+s1_s
            s2(ii)=s2(ii)+s2_s
         ENDDO
      ENDDO
!
      END SUBROUTINE sum_nf2_m
!
!------------------------------------------------------------------------------
!
      SUBROUTINE sum_nf23_nbcon(s1i_inl,s1i_out, &
                                s1)
!
      USE Zparam        , ONLY: nb_max
      USE Zvec_param    , ONLY: nf_inl,nf_out
      USE Znum_cell     , ONLY: istart_nb1,      &
                                ia_nb,           &
                                istart_nbcon_nf
      USE Zvec_index    , ONLY: nbcon_nf
!
      IMPLICIT NONE
!.....Input
      REAL(8) :: s1i_inl(nf_inl),s1i_out(nf_out)
!.....Output
      REAL(8) :: s1(nb_max)
!.....Local variables
      INTEGER :: i,k,nb
      INTEGER :: nf_number,len,istart1,i1
      INTEGER :: istart3,i3
!
      DO i=1,nb_max
         s1(i)=0.d0
      ENDDO
!
!.....Partial sum for inl
!
      nf_number=2
      istart3=istart_nbcon_nf(nf_number)
      istart1=istart_nb1(1,nf_number)
      len    =istart_nb1(2,nf_number)
      DO nb=1,len
         i1=istart1+nb
!DIR$ NOVECTOR
         DO i=ia_nb(i1),ia_nb(i1+1)-1
            i3=istart3+i
            k=nbcon_nf(i3)
            s1(k)=s1(k)-s1i_inl(i)
         ENDDO
      ENDDO
!
!.....Partial sum for out
!
      nf_number=3
      istart3=istart_nbcon_nf(nf_number)
      istart1=istart_nb1(1,nf_number)
      len    =istart_nb1(2,nf_number)
      DO nb=1,len
         i1=istart1+nb
!DIR$ NOVECTOR
         DO i=ia_nb(i1),ia_nb(i1+1)-1
            i3=istart3+i
            k=nbcon_nf(i3)
            s1(k)=s1(k)+s1i_out(i)
         ENDDO
      ENDDO
!
      END SUBROUTINE sum_nf23_nbcon
!
!------------------------------------------------------------------------------
!
      SUBROUTINE grad_press_ls_sum(eps,s1i_i,dx_non,dx_sym,s1,s2)
!
      USE Zzone         , ONLY: ncell_fluid
      USE Zparam        , ONLY: ndim
      USE Zvec_param    , ONLY: nf_non,nf_sym
      USE Znum_cell     , ONLY: istart_nf,istart_nb1,     &
                                ia_nb,icell_nb,right_nb_k
      USE Zbc_index    , ONLY: ngrad
      USE Zgrad_ls_c2d , ONLY: det_2
      USE Zgrad_ls_c3d , ONLY: lsindex,det_3
      USE Zvec_geo     , ONLY: dji_nf
!
      IMPLICIT NONE
!  
!.....Input
      REAL(8) :: eps
      REAL(8) :: s1i_i(nf_non,ndim)
      REAL(8) :: dx_non(nf_non,ndim)
      REAL(8) :: dx_sym(nf_sym,ndim)
!.....Output
      REAL(8) :: s1(ncell_fluid,ndim),s2(ncell_fluid,ndim)
!.....Local variables
      INTEGER :: i,k,nb
      INTEGER :: ii,ix
      INTEGER :: nf_number,istart,len,istart1,i1
      REAL(8) :: s1_s,s2_s
!.....Local arrays
      REAL(8) :: det(ncell_fluid)
!
      IF(ndim.eq.2) THEN
         DO i=1,ncell_fluid
            det(i)=det_2(i)
         ENDDO
      ELSE
         DO i=1,ncell_fluid
            det(i)=det_3(i)
         ENDDO
      ENDIF
!
      DO ix=1,ndim
!
!........Initial source values
!
         DO i=1,ncell_fluid
            s1(i,ix)=0.d0
            s2(i,ix)=0.d0
         ENDDO
!
!........Partial sum for non_k
!
         nf_number=-1
         len   =istart_nb1(2,nf_number)
         DO nb=1,len
            ii=icell_nb(nb)
            IF(ngrad(abs(ii)).le.0) cycle
            IF(lsindex(abs(ii)).eq.0) cycle
            IF(det(abs(ii)).eq.1.0) THEN
               s1_s=0.d0
!DIR$ NOVECTOR
               DO i=ia_nb(nb),ia_nb(nb+1)-1
                  k=right_nb_k(i)
                  s1_s=s1_s+s1i_i(k,ix)
               ENDDO
               s1(ii,ix)=s1(ii,ix)+s1_s
               s2(ii,ix)=0.d0
            ELSE
               s1_s=0.d0
               s2_s=0.d0
!DIR$ NOVECTOR
               DO i=ia_nb(nb),ia_nb(nb+1)-1
                  k=right_nb_k(i)
                  IF(dabs(dx_non(k,ix)).gt.eps)THEN
                     s1_s=s1_s+s1i_i(k,ix)/(dji_nf(k)**2)
                     s2_s=s2_s+1.d0
                  ENDIF
               ENDDO
               s1(ii,ix)=s1_s
               s2(ii,ix)=s2_s
            ENDIF
         ENDDO
!
!........Partial sum for non
!
         nf_number=0
         istart =istart_nf(1,nf_number)
         istart1=istart_nb1(1,nf_number)
         len    =istart_nb1(2,nf_number)
         DO nb=1,len
            i1=istart1+nb
            ii=icell_nb(i1)
            IF(ngrad(abs(ii)).le.0) cycle
            IF(lsindex(abs(ii)).eq.0) cycle
            IF(det(abs(ii)).eq.1.0) THEN
               s1_s=0.d0
!DIR$ NOVECTOR
               DO i=ia_nb(i1),ia_nb(i1+1)-1
                  s1_s=s1_s+s1i_i(i,ix)
               ENDDO
               s1(ii,ix)=s1(ii,ix)+s1_s
            ELSE
               s1_s=0.d0
               s2_s=0.d0
!DIR$ NOVECTOR
               DO i=ia_nb(i1),ia_nb(i1+1)-1
                  IF(dabs(dx_non(i,ix)).gt.eps)THEN
                     s1_s=s1_s+s1i_i(i,ix)/(dji_nf(i)**2)
                     s2_s=s2_s+1.d0
                  ENDIF
               ENDDO
               s1(ii,ix)=s1(ii,ix)+s1_s
               s2(ii,ix)=s2(ii,ix)+s2_s
            ENDIF
         ENDDO
!
!........Partial sum for sym
!
         nf_number=8
         istart =istart_nf(1,nf_number)
         istart1=istart_nb1(1,nf_number)
         len    =istart_nb1(2,nf_number)
         DO nb=1,len
            i1=istart1+nb
            ii=icell_nb(i1)
            IF(ngrad(abs(ii)).le.0) cycle
            IF(lsindex(abs(ii)).eq.0) cycle
            IF(det(abs(ii)).ne.1.0) THEN
               s2_s=0.d0
!DIR$ NOVECTOR
               DO i=ia_nb(i1),ia_nb(i1+1)-1
                  IF(dabs(dx_sym(i,ix)).gt.eps)THEN
                     s2_s=s2_s+1.d0
                  ENDIF
               ENDDO
               s2(ii,ix)=s2(ii,ix)+s2_s
            ENDIF
         ENDDO
!
      ENDDO
!
      END SUBROUTINE grad_press_ls_sum
!
!------------------------------------------------------------------------------
!
      SUBROUTINE grad_press_ls_c_2d_sum(s1_i,s2_i, &
                                        s1,s2,s3)
!
      USE Zzone        , ONLY: ncell_fluid
      USE Zparam       , ONLY: ndim
      USE Zvec_param   , ONLY: nf_non
      USE Znum_cell    , ONLY: istart_nb1,                &
                               ia_nb,icell_nb,right_nb_k
      USE Zbc_index    , ONLY: ngrad
      USE Zgrad_ls_c3d , ONLY: lsindex
!
      IMPLICIT NONE
!  
!.....Input
      REAL(8) :: s1_i(nf_non,ndim)
      REAL(8) :: s2_i(nf_non,ndim:ndim)
!.....Output
      REAL(8) :: s1(ncell_fluid),s2(ncell_fluid),s3(ncell_fluid)
!.....Local variables
      INTEGER :: i,k,nb
      INTEGER :: ii
      INTEGER :: nf_number,len,istart1,i1
      REAL(8) :: s1_s,s2_s,s3_s
!
!........Initial source values
!
      DO i=1,ncell_fluid
         s1(i)=0.d0
         s2(i)=0.d0
         s3(i)=0.d0
      ENDDO
!
!.....Partial sum for non_k
!
      nf_number=-1
      len   =istart_nb1(2,nf_number)
      DO nb=1,len
         ii=icell_nb(nb)
         IF(ngrad(abs(ii)).le.0) cycle
         IF(lsindex(abs(ii)).eq.0) cycle
            s1_s=0.d0
            s2_s=0.d0
            s3_s=0.d0
!DIR$ NOVECTOR
            DO i=ia_nb(nb),ia_nb(nb+1)-1
               k=right_nb_k(i)
               s1_s=s1_s+s1_i(k,1)
               s2_s=s2_s+s1_i(k,2)
               s3_s=s3_s+s2_i(k,2)
            ENDDO
               s1(ii)=s1(ii)+s1_s
               s2(ii)=s2(ii)+s2_s
               s3(ii)=s3(ii)+s3_s
      ENDDO
!
!.....Partial sum for non
!
      nf_number=0
      istart1=istart_nb1(1,nf_number)
      len    =istart_nb1(2,nf_number)
      DO nb=1,len
         i1=istart1+nb
         ii=icell_nb(i1)
         IF(ngrad(abs(ii)).le.0) cycle
         IF(lsindex(abs(ii)).eq.0) cycle
            s1_s=0.d0
            s2_s=0.d0
            s3_s=0.d0
!DIR$ NOVECTOR
            DO i=ia_nb(i1),ia_nb(i1+1)-1
               s1_s=s1_s+s1_i(i,1)
               s2_s=s2_s+s1_i(i,2)
               s3_s=s3_s+s2_i(i,2)
            ENDDO
            s1(ii)=s1(ii)+s1_s
            s2(ii)=s2(ii)+s2_s
            s3(ii)=s3(ii)+s3_s
      ENDDO
!
      END SUBROUTINE grad_press_ls_c_2d_sum
!
!------------------------------------------------------------------------------
!
      SUBROUTINE grad_press_ls_c_3d_sum(s1_i,s2_i,s3_i, &
                                        s1,s2,s3,s4,s5,s6)
!
      USE Zzone        , ONLY: ncell_fluid
      USE Zparam       , ONLY: ndim
      USE Zvec_param   , ONLY: nf_non
      USE Znum_cell    , ONLY: istart_nb1,               &
                               ia_nb,icell_nb,right_nb_k
      USE Zbc_index    , ONLY: ngrad
      USE Zgrad_ls_c3d , ONLY: lsindex
!
      IMPLICIT NONE
!  
!.....Input
      REAL(8) :: s1_i(nf_non,ndim)
      REAL(8) :: s2_i(nf_non,2:ndim)
      REAL(8) :: s3_i(nf_non,ndim:ndim)
!.....Output
      REAL(8) :: s1(ncell_fluid),s2(ncell_fluid),s3(ncell_fluid), &
                 s4(ncell_fluid),s5(ncell_fluid),s6(ncell_fluid)
!.....Local variables
      INTEGER :: i,k,nb
      INTEGER :: ii
      INTEGER :: nf_number,len,istart1,i1
      REAL(8) :: s1_s,s2_s,s3_s,s4_s,s5_s,s6_s
!
!........Initial source values
!
      DO i=1,ncell_fluid
         s1(i)=0.d0
         s2(i)=0.d0
         s3(i)=0.d0
         s4(i)=0.d0
         s5(i)=0.d0
         s6(i)=0.d0
      ENDDO
!
!.....Partial sum for non_k
!
      nf_number=-1
      len   =istart_nb1(2,nf_number)
      DO nb=1,len
         ii=icell_nb(nb)
         IF(ngrad(abs(ii)).le.0) cycle
         IF(lsindex(abs(ii)).eq.0) cycle
            s1_s=0.d0
            s2_s=0.d0
            s3_s=0.d0
            s4_s=0.d0
            s5_s=0.d0
            s6_s=0.d0
!DIR$ NOVECTOR
            DO i=ia_nb(nb),ia_nb(nb+1)-1
               k=right_nb_k(i)
               s1_s=s1_s+s1_i(k,1)
               s2_s=s2_s+s1_i(k,2)
               s3_s=s3_s+s1_i(k,3)
               s4_s=s4_s+s2_i(k,2)
               s5_s=s5_s+s2_i(k,3)
               s6_s=s6_s+s3_i(k,3)
            ENDDO
               s1(ii)=s1(ii)+s1_s
               s2(ii)=s2(ii)+s2_s
               s3(ii)=s3(ii)+s3_s
               s4(ii)=s4(ii)+s4_s
               s5(ii)=s5(ii)+s5_s
               s6(ii)=s6(ii)+s6_s
      ENDDO
!
!.....Partial sum for non
!
      nf_number=0
      istart1=istart_nb1(1,nf_number)
      len    =istart_nb1(2,nf_number)
      DO nb=1,len
         i1=istart1+nb
         ii=icell_nb(i1)
         IF(ngrad(abs(ii)).le.0) cycle
         IF(lsindex(abs(ii)).eq.0) cycle
            s1_s=0.d0
            s2_s=0.d0
            s3_s=0.d0
            s4_s=0.d0
            s5_s=0.d0
            s6_s=0.d0
!DIR$ NOVECTOR
            DO i=ia_nb(i1),ia_nb(i1+1)-1
               s1_s=s1_s+s1_i(i,1)
               s2_s=s2_s+s1_i(i,2)
               s3_s=s3_s+s1_i(i,3)
               s4_s=s4_s+s2_i(i,2)
               s5_s=s5_s+s2_i(i,3)
               s6_s=s6_s+s3_i(i,3)
            ENDDO
            s1(ii)=s1(ii)+s1_s
            s2(ii)=s2(ii)+s2_s
            s3(ii)=s3(ii)+s3_s
            s4(ii)=s4(ii)+s4_s
            s5(ii)=s5(ii)+s5_s
            s6(ii)=s6(ii)+s6_s
      ENDDO
!
      END SUBROUTINE grad_press_ls_c_3d_sum
