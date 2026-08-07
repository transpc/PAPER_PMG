include ../makefile.in

MYNAME = cupidMOD
MYNAME1 = s_cupidMOD
VER = 01
EXCNAME = $(MYNAME).x
LIBNAME = lib$(MYNAME).a
SRCNAME = $(MYNAME).for
BACKUPNAME = $(MYNAME1)_$(VER).tar

#----------------------------------------------------------------------
# extension type : .f90
FSRCS = cub_interp.f90\
	f_fric_face.f90\
	f_nd_face.f90\
	heat_partition_fluid.f90\
	heat_partition_2_bi_lo.f90\
	heat_partition_2_bi_lo_err.f90\
	IAT_bubsize_change.f90\
	IAT_calc.f90\
	IAT_source.f90\
	int_area.f90\
	int_drag.f90\
	int_htc.f90\
	int_htc_interp_temp.f90\
	int_htc_simple_model.f90\
	int_non_drag_coeff.f90\
	int_non_drag_lub.f90\
	int_non_drag_turb.f90\
	int_swap.f90\
	virtual_mass.f90\
	wall_drag.f90\
        wall_condensation_binary_diff.f90\
        wall_condensation_HMTA.f90\
        wall_condensation_model.f90\
        wall_condensation_RBLA.f90\
        wall_condensation_uchida.f90\
	wall_condensation_steam.f90\
	non_rv_wall_ht_2d.f90
FSRCV = int_htc_simple_topology.f90\
	int_htc_full_topology.f90\
	IAT_convection.f90\
	hibiki.f90\
	heat_partition.f90\
	heat_partition_porous.f90\
	radiation_model.f90\
#----------------------------------------------------------------------

OBJSS = $(FSRCS:.f90=.o)
OBJSV = $(FSRCV:.f90=.o)
OBJS = $(OBJSS) $(OBJSV)
$(OBJSS): FFLAGS := $(FOPTFLAGS) -module ../Modules
$(OBJSV): FFLAGS := $(FOPTFLAGV) -module ../Modules
FSRC = $(FSRC1)

BACKUPFILE = *.f90 makefile
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
	$(ARCH) $(ARCHFLAGS) $(LIBNAME) $(OBJS)
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

