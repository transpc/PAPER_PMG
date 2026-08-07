--------------------------------------------------------------------------------
CUPID needs pre-requisites as follows; 

      - Intel Fortran Compiler (higher than v15.0)
      - OpenMPI 3.1.0
      
CUPID libraries were compiled on such environment. In order to run parallel
computation with domain decomposition, CUPID uses METIS 5.0 as well.      

===========
Compilation
===========
Compiling CUPID requires "makefile" and "makefile.in" which is already given.
Options in "makefile.in" for compilation should be changed for the users' 
Hardware environment. The makefile.in has information as follows;
     
      - Fortran command (ifort, mpif90, mpiifort, etc.)
      - Optimazation option
      - Link to CUPID libraries
      - Link to METIS libraries
      
In order to compile CUPID, Execute the following command:

shell$ make

Executable file "cupid.x" is successfully generated.

Parallel compilations are also supported (although some versions of "make",
such as GNU make, will only use the first target listed on the command
line when executable parallel builds).  For example (assume GNU make):

shell$ make procs=XX (XX is number of processors)

===
Run
===
In order to run CUPID, It is recommended to copy cupid.x into run directory
where mandatory files are:

      - somaFlow.in (T/H input)
      - foamGrid.in (Grids)
      - tpfh2o      (Steamtable)
      - cupid.x     (CUPID executable)

For serial calculation;

shell$ ./cupid.x

For parallel calculation;

shell$ mpirun -np XX ./cupid.x

=================================
Post-processing (ONLY in Windows)
=================================
CUPID results can be visualized ONLY in Windows currently. 
Paraview (higher than 4.1.0) is recommended to visualize the CUPID results.
In the directory to be post-processed, mandatories are as follows;

      - Makeplot.exe
      - makeplot.in
      - foamGrid.in  (same as run directory)
      - somaPlot.viw (CUPID results)
      
Once execute "Makeplot.exe", new directory named "foam" will be generated. 
Inside that directory, open "run.foam" in the Paraview.