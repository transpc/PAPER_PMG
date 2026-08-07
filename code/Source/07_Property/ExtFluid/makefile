include ../../makefile.in

MYNAME = cupidPRO
MYNAME1 = s_cupidPRO
VER = 01
EXCNAME = $(MYNAME).x
LIBNAME = lib$(MYNAME).a
SRCNAME = $(MYNAME).for
BACKUPNAME = $(MYNAME1)_$(VER).tar

#----------------------------------------------------------------------
# extension type : .f90
FSRCS = std2x6.f90\
	std2x1.f90\
	std2x3.f90\
	nth2x1.f90\
	nth2x2.f90\
	nth2x3.f90\
	nth2x6f.f90\
	strpu.f90\
	strpu1.f90\
	strpu2.f90\
	strtp.f90\
	strtro.f90\
	strx.f90\
	strpx.f90\
	strtx.f90\
	strsat.f90\
	mixprops.f90\
	satprops.f90\
	setallprops.f90\
	setsngprops.f90\
	getsatatt.f90\
	getsatatp.f90\
	gettsat.f90\
	getsnglatpu.f90\
	getsnglatpt.f90\
	getindex.f90\
	setpointers.f90\
	setanchor.f90\
	moveloright.f90\
	movehileft.f90\
	movehiright.f90\
	moveloleft.f90\
	newstread.f90\
	loadcorners.f90\
	herm2dplus.f90\
	herm2d.f90\
	herm2dleft.f90\
	gfunct.f90\
	gfuncte.f90\
	gpfunct.f90\
	gpfuncte.f90\
	gpfunctline.f90\
	gtfunct.f90\
	gtfuncte.f90\
	gttfunct.f90\
	gttfuncte.f90\
	gttfunctline.f90\
	gppfunct.f90\
	gppfuncte.f90\
	gppfunctline.f90\
	gptfunct.f90\
	gptfuncte.f90\
	gptfunctline.f90\
	rholine.f90\
	kappaline.f90\
	betaline.f90\
	cpline.f90\
	beta2d.f90\
	cubic.f90\
	checkvalue.f90\
	clearprops.f90
FSRCV =
#----------------------------------------------------------------------

OBJSS = $(FSRCS:.f90=.o)
OBJSV = $(FSRCV:.f90=.o)
OBJS = $(OBJSS) $(OBJSV)
$(OBJSS): FFLAGS := $(FOPTFLAGS) -module ../../Modules
$(OBJSV): FFLAGS := $(FOPTFLAGV) -module ../../Modules

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

