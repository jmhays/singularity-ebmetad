# put here the dir in which you compile the baby
THISDIR = /Users/chicco/Programs/gronamd/md_meta_devel/utilities/catdcd_plumed
# put the arch as declared when you compile vmd plugins
MYARCH = MACOSXX86_64
# put the dir that contains the plugins 
MYPLUGINSDIR = /Users/chicco/Programs/gronamd/md_meta_devel/utilities/path_tools/HYBRID/plugins
# your tcl version and location
TCLLIB=-L/System/Library/Frameworks/ 
TCLLDFLAGS=-ltcl8.5
# compiler and flags
CXX =g++
CXXFLAGS =-g

# this should be all: hopefully you should not touch here below
INCDIR =-I${THISDIR} -I${THISDIR}/plumed_links  -I${MYPLUGINSDIR}/include -I${MYPLUGINSDIR}/molfile_plugin/src -I${MYPLUGINSDIR}/compile/lib_${MYARCH}/molfile

EXEC = catdcd_plumed 
cfiles = $(wildcard *.c) $(wildcard ./plumed_links/*.c) 
OBJ = $(cfiles:.c=.o) 

%.o: %.c
	${CXX} -DAMBER ${CXXFLAGS} ${INCDIR}  -c  -o $@ $*.c 

$(EXEC): $(OBJ) ${MYPLUGINSDIR}/compile/lib_${MYARCH}/molfile/libmolfile_plugin.a
	${CXX} ${CXXFLAGS} ${OBJ} -o $(EXEC) -L${MYPLUGINSDIR}/compile/lib_${MYARCH}/molfile -lmolfile_plugin ${TCLLIB} ${TCLLDFLAGS} -lm

clean:
	@rm -rf *.o catdcd plumed_links/*.o  

