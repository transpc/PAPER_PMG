!comdeck contrl
      COMMON/contrl_cupid/dthy,dtht,dtn,dt,timehy,timeht,errmax,tmass,tmasso, &
      emass,emasso,count,curtmi,curtmj,curtrs,prevnd,cpurem(5),stdtrn,  &
      gravcn,testda(20),aflag,isrestartenabled0,ismajoreditenabled1,    &
      isminoreditenabled2,isplotenabled3,iscompleterestart4,            &
      interactiveflag5,isimplicithttrnsfr6,istwostepflag7,              &
      stdystatetermflg8,integrationflag9,isathenaopt10,isscdapopt11,    &
      isplotsqzflg12,rwtrstpltcmprssflag13,transrotatflag14,            &
      rstblcknumber1531,succes,done,ncount,nstsp,nrepet,help,nany,skipt,&
      problemtype05,problemopt611,ncase,fail,uniti,unito,chngno(99),    &
      ishydrochng0,isheatcondchng1,isrknchng2,ispltrcrdchng3,           &
      isradiationchng4,isfissionprdctchng5,iscntrlsyschng6,pageno,      &
      isprntaccum0,isprntbrntrn1,isprntccfl2,isprntchfcal3,             &
      isprntconden4,isprntdittus5,isprnteqfinl6,isprntfwdrag7,          &
      isprntht2tdp9,isprnthtfinl12,isprnthtrc113,isprnthydro15,         &
      isprntistate17,isprntjchoke18,isprntjprop19,isprntnoncnd20,       &
      isprntphantj21,isprntphantv22,isprntpimplt23,isprntpintfc24,      &
      isprntprednb25,isprntpreseq26,isprntpstdnb27,isprntqfmove28,      &
      isprntsimplt29,isprntstacc0,isprntstate1,isprntstatep2,           &
      isprntsuboil3,isprntsysitr4,isprnttstate6,isprntvalve7,           &
      isprntvexplt8,isprntvfinl9,isprntvimplt10,isprntvlvela11,         &
      isprntvolvel12,isprnttrip13,isprntpower14,isprntvolume15,         &
      isprntjunction16,isprntheatstr17,isprntradht18,isprntreflood19,   &
      isprntfsnprdtr20,isprntcontrol21,isprntinput22,isprntmiedit23,    &
      islevelmodel0,islevelcrossing1,issnapon,safe2
      !gamma rktpow3d has been removed in MARS. The variable belongs to 
      !PARCS code.
      REAL(8) dthy,dtht,dtn,dt,timehy,timeht,errmax,tmass,tmasso,emass, &
      emasso,count,curtmi,curtmj,curtrs,prevnd,cpurem,stdtrn,gravcn,    &
      testda,safe2
      INTEGER succes,done,ncount,nstsp,nrepet,help,nany,ncase,pageno,   &
      cpurei(2,5)
      LOGICAL aflag,skipt,fail,uniti,unito,chngno,nmechk,issnapon
      !gamma
      ! replacement variables for iextra
      LOGICAL ishydrochng0,isheatcondchng1,isrknchng2,ispltrcrdchng3,   &
      isradiationchng4,isfissionprdctchng5,iscntrlsyschng6 
      !
      ! replacement variables for imdctl
      LOGICAL islevelmodel0,islevelcrossing1 
      !
      ! replacement variables for iroute
      INTEGER problemtype05,problemopt611 
      !
      ! replacement variables for print
      LOGICAL isrestartenabled0,ismajoreditenabled1,isminoreditenabled2,&
      isplotenabled3,iscompleterestart4,interactiveflag5,               &
      isimplicithttrnsfr6,istwostepflag7,stdystatetermflg8,             &
      integrationflag9,isathenaopt10,isscdapopt11,isplotsqzflg12,       &
      rwtrstpltcmprssflag13,transrotatflag14
      !
      ! replacement variables for ihlppr
      LOGICAL isprntaccum0,isprntbrntrn1,isprntccfl2,isprntchfcal3,     &
      isprntconden4,isprntdittus5,isprnteqfinl6,isprntfwdrag7,          &
      isprntht2tdp9,isprnthtfinl12,isprnthtrc113,isprnthydro15,         &
      isprntistate17,isprntjchoke18,isprntjprop19
      !
      LOGICAL isprntnoncnd20,isprntphantj21,isprntphantv22,             &
      isprntpimplt23,isprntpintfc24,isprntprednb25,isprntpreseq26,      &
      isprntpstdnb27,isprntqfmove28,isprntsimplt29
      !
      LOGICAL isprntstacc0,isprntstate1,isprntstatep2,isprntsuboil3,    &
      isprntsysitr4,isprnttstate6,isprntvalve7,isprntvexplt8,           &
      isprntvfinl9,isprntvimplt10,isprntvlvela11,isprntvolvel12,        &
      isprnttrip13,isprntpower14,isprntvolume15 
      !
      LOGICAL isprntjunction16,isprntheatstr17,isprntradht18,           &
      isprntreflood19,isprntfsnprdtr20,isprntcontrol21,isprntinput22,   &
      isprntmiedit23
      !
      INTEGER rstblcknumber1531 
      !gammaend
      EQUIVALENCE(safe2,nmechk),(cpurem(1),cpurei(1,1)) 
!dab+
!  Global variables for storing restart edit summary information.       
      INTEGER mxrst 
      PARAMETER (mxrst=500) 
      INTEGER nrstblk,irstno(mxrst),iblkno(mxrst) 
      REAL(8) rtime(mxrst) 
      COMMON/rstsum/rtime,nrstblk,irstno,iblkno 
!dab-
!
!***********************************************************************
!
! Data dictionary for local variables
!
! Number of local variables =  116
!
! i=integer r=real l=logical c=character
!***********************************************************************
! Type    Name                              Definition
!-----------------------------------------------------------------------
!  l aflag                 = Set true in dtstep when heat
!                           structures are to be advanced.
!                           checked in htadv
!  l chngno(90)            = true if option # is set
!  r count                 = Real value of the integer ncount
!  i cpurei(5)             = Integer values (5) of the real cpurem
!                           (1,2,3) not used
!                           (4) word 4 from 105 card
!                            = ncount1, ncount for starting
!                               diagnostic edits
!                            = -1, write to dumpfil1 for
!                               ncount2 = cpurei(5), stop
!                            = -2, write to dumpfil1 for
!                               ncount2 = cpurei(5),
!                               redo timestep, write to
!                               dumpfil2, stop
!                            = -3, set in dtstep for -2 case.
!                           (5) word 5 from 105 card
!                            = ncount2, ncount for
!                               terminating diagnostic edits
!  r cpurem(5)             = Contains cpu remaining times values
!                            and advancement counts to start/stop
!                            diagnostic edit from input.
!  r curtmi                = Time for next minor edit.                  
!  r curtmj                = Time for next major edit.                  
!  r curtrs                = Time for next restart edit.                
!  i done                  = Integer flag indicating state of
!                            transient.
!                            0   = advancements are to continue,
!                            !=0 = advancements are to be
!                             terminated.
!  r dt                    = Current time step.
!  r dtht                  = Heat structure time step, currently
!                            the same as dthy.
!  r dthy                  = Hydrodynamic time step requested by
!                            user.
!  r dtn                   = Time step limit due to material
!                            transport stability limit
!                           (Courant limit).
!  r emass                 = Current mass error.
!  r emasso                = Old mass error.
!  r errmax                = Error estimate used in time step
!                            control.
!  l fail                  = Set to true if error encountered.
!  r gravcn                = Gravitational constant (which may be
!                            set by input).
!  i help                  = Used to control debug editing and
!                            termination.
!  l integrationflag9      = on-line selection of time integration
!                            flag
!  l interactiveflag5      = interactive flag.
!  l isathenaopt10         = athena option.
!  l iscntrlsyschng6       = control system change.
!  l iscompleterestart4    = complete restart
!  l isfissionprdctchng5   = fission product change.
!  l isheatcondchng1       = heat conduction change.
!  l ishydrochng0          = hydrodynamic change.
!  l isimplicithttrnsfr6   = implicit heat transfer
!  l islevelcrossing1      = level crossing activated.
!  l islevelmodel0         = level model activated.
!  l ismajoreditenabled1   = major edit enabled
!  l isminoreditenabled2   = minor edit enabled
!  l isplotenabled3        = plot enabled
!  l isplotsqzflg12        = plot record squoz flag.
!  l ispltrcrdchng3        = plot record change.
!  l isprntaccum0          = activate accum print block
!  l isprntbrntrn1         = activate brntrn print block
!  l isprntccfl2           = activate ccfl print block
!  l isprntchfcal3         = activate chfcal print block
!  l isprntconden4         = activate conden print block
!  l isprntcontrol21       = activate control print block
!  l isprntdittus5         = activate dittus print block
!  l isprnteqfinl6         = activate eqfinl print block
!  l isprntfsnprdtr20      = activate fsnprdtr print block
!  l isprntfwdrag7         = activate fwdrag print block
!  l isprntheatstr17       = activate heatstr print block
!  l isprntht2tdp9         = activate ht2tdp print block
!  l isprnthtfinl12        = activate htfinl print block
!  l isprnthtrc113         = activate htrc print block
!  l isprnthydro15         = activate hydro print block
!  l isprntinput22         = activate input print block
!  l isprntistate17        = activate istate print block
!  l isprntjchoke18        = activate jchoke print block
!  l isprntjprop19         = activate jprop print block
!  l isprntjunction16      = activate junction print block
!  l isprntmiedit23        = activate miedit print block
!  l isprntnoncnd20        = activate noncnd print block
!  l isprntphantj21        = activate phantj print block
!  l isprntphantv22        = activate phantv print block
!  l isprntpimplt23        = activate phantv pimplt block
!  l isprntpintfc24        = activate phantv pintfc block
!  l isprntpower14         = activate power print block
!  l isprntprednb25        = activate prednb print block
!  l isprntpreseq26        = activate preseq print block
!  l isprntpstdnb27        = activate pstdnb print block
!  l isprntqfmove28        = activate qfmove print block
!  l isprntradht18         = activate radht print block
!  l isprntreflood19       = activate reflood print block
!  l isprntsimplt29        = activate simplt print block
!  l isprntstacc0          = activate stacc print block
!  l isprntstate1          = activate state print block
!  l isprntstatep2         = activate statep print block
!  l isprntsuboil3         = activate suboil print block
!  l isprntsysitr4         = activate sysitr print block
!  l isprnttrip13          = activate trip print block
!  l isprnttstate6         = activate tstate print block
!  l isprntvalve7          = activate valve print block
!  l isprntvexplt8         = activate vexplt print block
!  l isprntvfinl9          = activate vfinl print block
!  l isprntvimplt10        = activate vimplt print block
!  l isprntvlvela11        = activate vlvela print block
!  l isprntvolume15        = activate volume print block
!  l isprntvolvel12        = activate volvel print block
!  l isradiationchng4      = radiation change.
!  l isrestartenabled0     = restart enable bit.
!  l isrknchng2            = reactor kinetics change.
!  l isscdapopt11          = scdap is active
!  l issnapon              = flag set to true if the snap socket is     
!                            active                                     
!  l istwostepflag7        = two step method
!  i nany                  = Number of mass error messages
!                            remaining to be issued.
!  i ncase                 = Case number
!                           Initially = -1 in blkdta
!                           0 if last case of series
!  i ncount                = Count of number of advancements,
!                            successful or otherwise.
!  l nmechk                = true if no mass error check requested
!                            on time step card
!  i nrepet                = Number of hydrodynamic advancements
!                            at current dt to finish a requested
!                            time step of dthy.
!  i nstsp                 = Number of standard advancements.
!  i pageno                = Page number.
!  r prevnd                = Previous time step card end time.          
!  i problemopt611         = Problem option (see ityp2 in inputd).
!                           1 = stdy-st (1 and 2 are used for new
!                                 and restart)
!                           2 = transnt
!                           3 = binary (3 and 4 are used for strip
!                           4 = fmtout
!  i problemtype05         = Problem type (see ityp1 in inputd).
!                            1 = new
!                            2 = restart
!                            3 = relap5 internal plots (not used)
!                            5 = strip
!                            6 = cmpcom (compare dump records)
!  i rstblcknumber1531     = restart block number.
!  l rwtrstpltcmprssflag13 =
!  r safe2                 = Last word in common block. Its address is
!                            obtained using locf and is used with the
!                            address of the first word, dthy, to
!                            compute the length of the common block.
!                            The length is used when the common block
!                            is written out to a file.
!                            Make sure it is on an 8-byte boundary.
!  l skipt                 = Logical flag used to control first
!                            entry processing in
!                            subroutine dtstep.
!  r stdtrn                = Steady state - transient flag
!  l stdystatetermflg8     = steady state termination
!  i succes                = Integer flag indicating success of
!                            the advancement.
!                            0 = no need to repeat advancement
!                                 with reduced time step
!                            1 = excessive truncation error
!                            2 = water property error
!                            3 = non-diagonal matrix
!                            4 = metal appears
!  r testda(20)            = Data array for minor edits and
!                            plotting during debug.
!                           Temporary coding is used to load data
!                            in this array and scnreq accepts
!                            testda as a legal request.
!  r timeht                = Problem time for heat structure
!                            advancements.
!  r timehy                = Problem time for hydrodynamic
!                            advancements.
!  r tmass                 = Total mass of water in system
!                            currently.
!  r tmasso                = Total mass of water in system at
!                            time = 0.0.
!  l transrotatflag14      = transient rotation flag
!                            0.0 = steady state
!                            1.0 = transient
!  l uniti                 = If true, SI units on input.
!  l unito                 = If true, SI units on output.
!***********************************************************************
!