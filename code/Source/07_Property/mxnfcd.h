!comdeck mxnfcd                                                         
!                                                                       
      INTEGER mxnfld 
      PARAMETER (mxnfld=17) 
!                                                                       
      COMMON/mxnfbv_cupid/wmoles(mxnfld) 
      REAL(8) wmoles 
!                                                                       
      COMMON/mxnfcv_cupid/tpfnam(mxnfld),tpfncl(mxnfld),tpfnin(mxnfld),fsymbl(&
      mxnfld)                                                           
      CHARACTER tpfnam * 256,tpfncl * 40,tpfnin * 40,fsymbl * 8 
