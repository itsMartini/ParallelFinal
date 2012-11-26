CC = mpif90
OPT = -O3
LIB = -llapack
MODS = variables_mod.mod ode_solver_mod.mod solutions_mod.mod w_mod.mod output_mod.mod
MOD_OBJS = variables_mod.o ode_solver_mod.o solutions_mod.o w_mod.o output_mod.o
F90 = variables_mod.f90 ode_solver_mod.f90 solutions_mod.f90 w_mod.f90 output_mod.f90 driver.f90 
OBJS = driver.o

all: $(MOD_OBJS) $(OBJS)
	$(CC) $(OPT) $(MOD_OBJS) $(OBJS) $(LIB) -o driver_exe

variables_mod.o: variables_mod.f90
	$(CC) $(OPT) -c variables_mod.f90

solutions_mod.o: solutions_mod.f90
	$(CC) $(OPT) -c solutions_mod.f90

w_mod.o: w_mod.f90
	$(CC) $(OPT) -c w_mod.f90

output_mod.o: output_mod.f90
	$(CC) $(OPT) -c output_mod.f90

ode_solver_mod.o: ode_solver_mod.f90
	$(CC) $(OPT) -c ode_solver_mod.f90

driver.o: driver.f90
	$(CC) $(OPT) -c driver.f90

clean:
	rm *.mod *.o *_exe *~