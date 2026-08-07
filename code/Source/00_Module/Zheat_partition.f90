      MODULE Zheat_partition
! 
      IMPLICIT NONE
      SAVE
!
      REAL(8) q1,qq,qe,qcl,qcg,d_depart,ndensity,bfreq,a_two,dlo &
              ,twait,hconvl,hconvg,kfactor,utaul
      REAL(8) q1ary(3),qqary(3),qeary(3),qclary(3),qcgary(3),weightary(3),&
              d_departary(3),ndensityary(3),bfreqary(3),a_twoary(3),dloary(3),errary(3),deltatsary(3)
!
      END MODULE Zheat_partition