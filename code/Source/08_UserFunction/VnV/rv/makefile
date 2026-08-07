include ../../../makefile.in

MYNAME = cupidUSER
MYNAME1 = s_cupidUSER
VER = 01
EXCNAME = $(MYNAME).x
LIBNAME = lib$(MYNAME).a
SRCNAME = $(MYNAME).for
BACKUPNAME = $(MYNAME1)_$(VER).tar

# LFLAGS = -L/cm/shared/apps/metis/METIS -lmetis 

FFLAGS = $(FOPTFLAGS) -module ../../../Modules
# CFLAGS = -I$(COPTFLAGS) \-L$(INCDIR)
CFLAGS =
#----------------------------------------------------------------------
# extension type : .f90
FSRC1 = fs_1d_out.f90\
	rbht_1d_out.f90\
	apr1400_lbloca_user.f90\
	apr1400_lbloca_out.f90\
	opr1000_rv_user.f90\
	geodata_user_opr1000_1x1.f90\
	geodata_user_opr1000_rod.f90\
	rv_wall_ht_2d_rod.f90\
	rv_wallht_nucl_boiling_HTC_rod.f90\
	rv_wallht_single_phase_HTC_rod.f90\
	rv_wallht_subcooled_boiling_rod.f90\
	rv_wallht_condensation_HTC_rod.f90\
	rv_wallht_trans_film_boiling_HTC_rod.f90\
	rv_wallht_trans_film_reflood_HTC_rod.f90\
	udfn_subchannel_info.f90\
	udfn_subchannel_chf.f90\
	apr1400_lbloca_ctrl.f90\
	udfn_mom_source_evvd_rod.f90\
	geodata_user_apr1400_1x1.f90
#----------------------------------------------------------------------

FSRC = $(FSRC1)

SRC  = $(FSRC)

OBJS = $(FSRC1:.f90=.o)

BACKUPFILE = *.f90 *.h makefile
.SUFFIXES: .f90 .o
.F.o :
	$(FC) $(FFLAGS) -c $<
.f.o :
	$(FC) $(FFLAGS) -c $<
.f90.o :
	$(FC) $(FFLAGS) -c $<

default : mk_lib

mk_exe: $(OBJS)
	$(FC) -o $(EXCNAME) $(OBJS) $(LFLAGS)

mk_lib: $(OBJS)

clean:
#	rm -f *~
#	rm -f *.o 
#	rm -f *.mod
#	rm -f *__genmod.*
#	rm -f $(EXCNAME)
#	rm -f $(LIBNAME)
#	rm -f $(SRCNAME)
	
new : clean default

src :
	cat $(SRC) > $(SRCNAME)
tag :
	ctags $(SRC)
backup :
	tar cvf $(BACKUPNAME) $(BACKUPFILE)
