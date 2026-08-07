mv makefile.lib1 makefile
mv ./06_MPI/makefile.lib ./06_MPI/makefile
make
find 00_Module -name '*.f90' -exec rm -f {} \;
find 00_Module -name '*.o' -exec rm -f {} \;
find 00_Module -name 'make*' -exec rm -f {} \;
find 03_Model -name '*.f90' -exec rm -f {} \;
find 03_Model -name '*.o' -exec rm -f {} \;
find 03_Model -name 'make*' -exec rm -f {} \;
find 04_Turbulence -name '*.f90' -exec rm -f {} \;
find 04_Turbulence -name '*.o' -exec rm -f {} \;
find 04_Turbulence -name 'make*' -exec rm -f {} \;
find 05_Solver -name '*.f90' -exec rm -f {} \;
find 05_Solver -name '*.o' -exec rm -f {} \;
find 05_Solver -name 'make*' -exec rm -f {} \;
find 06_MPI -name 'all*.f90' -exec rm -f {} \;
find 06_MPI -name 'bca*.f90' -exec rm -f {} \;
find 06_MPI -name 'com*.f90' -exec rm -f {} \;
find 06_MPI -name '*.o' -exec rm -f {} \;
find 07_Property -name '*.f90' -exec rm -f {} \;
find 07_Property -name '*.o' -exec rm -f {} \;
find 07_Property -name 'make*' -exec rm -f {} \;
find 10_LinkToMARS -name '*.f90' -exec rm -f {} \;
find 10_LinkToMARS -name '*.o' -exec rm -f {} \;
find 10_LinkToMARS -name 'make*' -exec rm -f {} \;
find 13_Vectorize -name '*.f90' -exec rm -f {} \;
find 13_Vectorize -name '*.o' -exec rm -f {} \;
find 13_Vectorize -name 'make*' -exec rm -f {} \;
find Closed -name '*.f90' -exec rm -f {} \;
find Closed -name '*.o' -exec rm -f {} \;
find Closed -name 'make*' -exec rm -f {} \;
rm -rf 11_LinkToMASTER
rm compile_cupid
mv makefile.in.lib makefile.in
mv makefile.lib2 makefile
mv ./06_MPI/makefile.exe ./06_MPI/makefile
