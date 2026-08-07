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
FSRC1 = tdiag.f90\
	udfn_1d_energy.f90\
	udfn_1d_momentum.f90\
	udfn_convert_t2e_liq.f90\
	udfn_heat_tube.f90\
	udfn_prop_in.f90\
	udfn_prop_liquid.f90\
	udfn_sg_allocate.f90\
	udfn_sg_heat_partition.f90\
	udfn_sg_heat_source.f90\
	udfn_sg_htc.f90\
	udfn_sg_input.f90\
	udfn_sg_mom_wall.f90\
	udfn_sg_outlet_property.f90\
	udfn_sg_post.f90\
	udfn_sg_pre.f90\
	udfn_sg_primary.f90\
	udfn_sg_primary_ini.f90\
	udfn_sg_thcond_lw.f90\
	udfn_sg_viscos_lw.f90\
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
