include ../makefile.in

MYNAME = cupidUSER
MYNAME1 = s_cupidUSER
VER = 01
EXCNAME = $(MYNAME).x
LIBNAME = lib$(MYNAME).a
SRCNAME = $(MYNAME).for
BACKUPNAME = $(MYNAME1)_$(VER).tar

#----------------------------------------------------------------------
# extension type : .f90
FSRCS = initialize_least_square_option.f90\
	initialize_specific_variables.f90\
	initialize_specific_variables_solid.f90\
	initialize_time_variant_inlet.f90\
	initialize_topology.f90\
	initialize_topology_criteria.f90\
	udfn_check_Hik.f90\
	udfn_flashing_hif.f90\
	udfn_mass_diff_coeff.f90\
	udfn_mat_prop.f90\
	udfn_mom_drag_i.f90\
	udfn_mom_source.f90\
	udfn_mom_wall.f90\
	udfn_outlet_property.f90\
	udfn_porous_hyd.f90\
	udfn_porous_property.f90\
	udfn_porous_user.f90\
	udfn_set_qvol_porous.f90\
	udfn_subchan_scalar_conv_mv.f90\
	udfn_turb_disp_model.f90\
	udfn_turb_zero.f90\
	udfn_vel_bc.f90\
	udfn_wall_drag_model.f90\
	udfn_weight.f90\
	udfn_mom_film_shear.f90\
	udfn_mom_press_source.f90\
	udfn_grid_user.f90\
	udfn_model_overwrite.f90\
	udfn_horizontal_outlet.f90
FSRCV = udfn_calc_HTC_int_i.f90\
	udfn_update_scalar.f90\
	udfn_erg_diff.f90\
	udfn_wall_lub_model.f90\
	udfn_heat_wallHTC_porous.f90\
#----------------------------------------------------------------------

OBJSS = $(FSRCS:.f90=.o)
OBJSV = $(FSRCV:.f90=.o)
OBJS = $(OBJSS) $(OBJSV)
$(OBJSS): FFLAGS := $(FOPTFLAGS) -module ../Modules
$(OBJSV): FFLAGS := $(FOPTFLAGV) -module ../Modules

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
