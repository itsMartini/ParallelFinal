Parallel Scientific Computing - Final Project - Kyle Geyser, Dylan Denning - README
------------------------------------------------------------------------------------
All of the results we obtained can be found in the directory results/. To produce
any results, first edit the file driver.f90 and uncomment any commented out lines
in the appropriate block (each block has a header describing its purpose).

Once the desired blocks have been uncommented from driver.f90, return to the command
line and type
     
	make

to compile all files into an executable named driver_exe. To run the program with one 
core, type

	./driver_exe

To run the program on a Sayers Lab machine, make sure you have an appropriate 
machine file in the current directory and type

	mpiexec -machinefile name_of_machine_file -np xx driver_exe

where xx is the number of cores you wish to use

To run the program on Mio, adjust the file my_project.pbs according to how many
nodes and cores per node you would like to use, then type

      qsub my_project.pbs

Note: for parallel performance testing, the parameters N and m are read from a 
file named config.txt (in that order). Modify this text file as necessary to 
change N and m.

To generate the plots shown in the report, run the python script plots.py either
by typing

      python plots.py

or type
      
      chmod u+x plots.py
then
      ./plots.py 

Note: mencoder must be installed to generate the movie (comment out appropriate lines
to ignore movie creation)

To clean the directory of all mod files, executables, auto-generated backup files,
etc., type

      make clean
------------------------------------------------------------------------------------
Detailed information on the project and our implementation of various algorithms can 
be found in the file FinalReport.pdf.

As mentioned above, all results can be found in the directory results/. Specifically,
results/ contains three tables created for two different examples (described in detail
in the report). results/mio contains the file runtimes.txt which lists N, m, number of
cores, and runtimes in a tabular format for Mio. results/sayers contains the file 
runtimes.txt which lists N, m, number of cores, and runtimes in a tabular format for
Sayers Lab. reults/plots/ contains all of the parallel performance plots as well as 
the surface plots we generated for analyzing our data. results/movies/ contains an mpg
file of an animation. results/tables/ is an empty directory used for storing formatted
output for LaTeX code.

Output of tables and runtimes upon running of the code will appear in results/.
Runtime data is appended to results/runtimes.txt.
