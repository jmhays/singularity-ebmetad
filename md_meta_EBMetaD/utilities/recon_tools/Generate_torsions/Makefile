.PHONY: all clean
all: GEN_TORSIONS

GEN_TORSIONS: gen_torsions.o generate_torsions.o
	gfortran generate_torsions.o gen_torsions.o -lstdc++ -o GEN_TORSIONS

gen_torsions.o: gen_torsions.cpp gen_torsions.h residue_library.h
	g++ -c -g gen_torsions.cpp -o gen_torsions.o

generate_torsions.o: generate_torsions.f90
	gfortran -c -g generate_torsions.f90 -o generate_torsions.o

clean:
	rm *.o GEN_TORSIONS

