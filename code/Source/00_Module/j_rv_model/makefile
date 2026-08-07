include ../../makefile.in

MYNAME = RvMODULE
MYNAME1 = s_RvMODULE
VER = 01
EXCNAME = $(MYNAME).x
LIBNAME = lib$(MYNAME).a
SRCNAME = $(MYNAME).for
BACKUPNAME = $(MYNAME1)_$(VER).tar

FFLAGS = $(FOPTFLAGS) -module ../../Modules 
CFLAGS = 
#----------------------------------------------------------------------
# extension type : .f90
FSRC1 = Zrv_flowmap.f90\
	Zrv_htc.f90\
	Zrv_hts_1d.f90\
	Zrv_hts_2d.f90\
	Zrv_int_friction.f90\
	Zrv_model.f90\
	Zrv_mpi.f90\
	Zrv_ncell.f90\
	Zrv_weight.f90\
	Zwall_HTC.f90\
	Zrv_ihtc_models.f90\
	Zrv_ihtc_pkg.f90\
	Zrv_wall_fric.f90\
	Zrv_choke.f90\
	Zrv_gap_cond.f90\
	Zrv_subchan.f90\
	Zrad_comp.f90
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
	$(FC) -o $(EXCNAME) $(OBJS) $(CFLAGS) $(LFLAGS)

mk_lib: $(OBJS)
	$(ARCH) $(ARCHFLAGS) $(LIBNAME) $(OBJS)
#
Zrv_ihtc_models.o: Zrv_ihtc_models.f90 Zrv_htc.o Zwall_HTC.o
Zrv_ihtc_pkg.o: Zrv_ihtc_pkg.f90 Zrv_htc.o Zrv_ihtc_models.o Zwall_HTC.o
#
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
