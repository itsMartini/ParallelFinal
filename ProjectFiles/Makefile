CC = gfortran
OPT = -O3
MODS = variables_mod.mod
MOD_OBJS = variables_mod.o
F90 = variables_mod.f90

all: $(MOD_OBJS)

variables_mod.o: variables_mod.f90
	$(CC) $(OPT) -c variables_mod.f90

clean:
	rm *.mod *.o *_exe *~