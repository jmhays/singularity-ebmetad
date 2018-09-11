# Note to self:
# To patch must patch everything in $plumedir/common_files and
# $plumedir/recon_plug/reconplug.a + $plumedir/recon_plug/src/recon_cbind.h

.PHONY: all clean
all: TRAJ_GREP 

include ./plumed.inc

TRAJ_GREP: traj_grep.o search.o $(PLUMED_OBJECTS) 
	gfortran $(PLUMED_OBJECTS) search.o traj_grep.o -llapack -lstdc++ -o TRAJ_GREP 

%.o: %.c
	gcc -c $< -g -DRECON_DRIVER
%.o: %.cpp
	g++ -c $< -g -DRECON_DRIVER

search.o: search.C
	g++ -c -g -DRECON_DRIVER search.C -o search.o 

traj_grep.o: traj_grep.f90
	gfortran -c -g traj_grep.f90 -o traj_grep.o

clean:	
	rm *.o TRAJ_GREP 
