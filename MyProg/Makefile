#
#   Makefile
#
FC = gfortran
FFLAGS = -Og -std=f95 --pedantic -ftrapv -fcheck=all -g -fbacktrace -ffpe-trap=zero,invalid,overflow -fcheck=bound   
LIBS = -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib -lblas -llapack 
#LIBS = -L/opt/OpenBLAS/lib -lopenblas

MODS   = diagonalization.o
OBJS   = main.o
#
all:    $(MODS) $(OBJS)
	$(FC) $(FFLAGS) -o main.exe $(OBJS) $(MODS) $(LIBS)
#
%.o: %.f
	$(FC) $(FFLAGS) -c $*.f
%.o: %.f90
	$(FC) $(FFLAGS) -c $*.f90
#
clean:
	rm -fr $(MODS) $(OBJS) *.exe *.mod
