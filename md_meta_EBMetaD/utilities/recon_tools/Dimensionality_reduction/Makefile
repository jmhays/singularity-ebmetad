# Note to self:
# To patch must patch everything in $plumedir/common_files and
# $plumedir/recon_plug/reconplug.a + $plumedir/recon_plug/src/recon_cbind.h

.PHONY: all clean
all: DIM_REDUCE 

RECON_OBJECTS= recon_basins.o \
               recon_utils.o \
               recon_cbind.o \
               recon_bespoke.o 

PLUMED_OBJECTS= \
        biasexchange.o \
        hills.o \
        metadyn.o \
        ptmetad.o \
        read_restraint.o \
        restraint_alfabeta.o \
        restraint_angle.o \
        restraint.o \
        restraint_coord.o \
        restraint_dihcor.o \
        restraint_dipole.o \
        restraint_dist.o \
        restraint_elstpot.o \
        restraint_hbonds.o \
        restraint_mindist.o \
        restraint_position.o \
        restraint_puckering.o \
        restraint_rgyr.o \
        restraint_rmsdtor.o \
        restraint_spath.o \
        restraint_torsion.o \
        restraint_waterbridge.o \
        restraint_zpath.o \
        restraint_helix.o \
        restraint_energy.o \
        restraint_alpharmsd.o \
        restraint_parabetarmsd.o \
        restraint_antibetarmsd.o \
        restraint_cmap.o \
        restraint_bespoke.o \
        testderivatives.o


DIM_REDUCE: reduce_dim.o dim_reduce.o $(RECON_OBJECTS) $(PLUMED_OBJECTS)  
	gfortran $(RECON_OBJECTS) $(PLUMED_OBJECTS) reduce_dim.o dim_reduce.o -llapack -lstdc++ -o DIM_REDUCE 

%.o: %.c
	gcc -c $< -g -DRECON_DRIVER
%.o: %.cpp
	g++ -c $< -g -DDEBUG -DRECON_DRIVER

dim_reduce.o: dim_reduce.C 
	g++ -c -g -DDEBUG -DRECON_DRIVER dim_reduce.C -o dim_reduce.o 

reduce_dim.o: reduce_dim.f90
	gfortran -c -g reduce_dim.f90 -o reduce_dim.o

clean:	
	rm *.o DIM_REDUCE 
