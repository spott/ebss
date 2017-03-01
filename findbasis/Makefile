include ../makefile.include

# To generate .clang_complete:
# make make PCC='~/.vim/bundle/clang_complete/bin/cc_args.py /opt/local/petsc-3.3-gcc-debug/bin/mpicxx' -B <target>

#OMPI_CXX=g++

LDFLAGS+= -lboost_mpi-mt -lboost_serialization-mt
#-lboost_mpi_mt
#OMPI_CXX=clang++-HEAD
TESTS=bspline/bspline-test.cpp
SOURCES=findbasis.cpp
OBJECTS=$(SOURCES:.cpp=.o)
TOBJECTS=$(TESTS:.cpp=.o)
EXECUTABLE=findbasis
TEXECUTABLE=bspline/bspline-test
DEFAULT=all


all: $(SOURCES) $(EXECUTABLE)

${TEXECUTABLE}: ${TOBJECTS}  chkopts style chk_git
	-${CLINKER}  -Wl,-keep_dwarf_unwind -o ${TEXECUTABLE} ${TOBJECTS} ${LDFLAGS} ${PETSC_VEC_LIB} ${SLEPC_LIB}

${EXECUTABLE}: ${OBJECTS}  chkopts style chk_git
	-${CLINKER}  -Wl,-keep_dwarf_unwind -o ${BINDIR}/${EXECUTABLE} ${OBJECTS} ${LDFLAGS} ${PETSC_VEC_LIB} ${SLEPC_LIB}

debug-ground: CPP_FLAGS+=-DDEBUG -DDEBUGEND -DDEBUGFINAL -DDEBUGGROUND1 -DDEBUGGROUND
debug-ground:  ${OBJECTS}  chkopts style chk_git
	-${CLINKER}  -Wl,-keep_dwarf_unwind -o ${BINDIR}/${EXECUTABLE} ${OBJECTS} ${LDFLAGS} ${PETSC_VEC_LIB} ${SLEPC_LIB}

debug-bound: CPP_FLAGS+=-DDEBUG -DDEBUGEND -DDEBUGFINAL -DDEBUGBOUND
debug-bound:  ${OBJECTS}  chkopts style chk_git
	-${CLINKER}  -Wl,-keep_dwarf_unwind -o ${BINDIR}/${EXECUTABLE} ${OBJECTS} ${LDFLAGS} ${PETSC_VEC_LIB} ${SLEPC_LIB}

debug-nodes: CPP_FLAGS+=-DDEBUG -DDEBUGEND -DDEBUGFINAL -DDEBUGNODES
debug-nodes:  ${OBJECTS}  chkopts style chk_git
	-${CLINKER}  -Wl,-keep_dwarf_unwind -o ${BINDIR}/${EXECUTABLE} ${OBJECTS} ${LDFLAGS} ${PETSC_VEC_LIB} ${SLEPC_LIB}

debug-excited: CPP_FLAGS+=-DDEBUG -DDEBUGEND -DDEBUGFINAL -DDEBUGEXCITED
debug-excited:  ${OBJECTS}  chkopts style chk_git
	-${CLINKER}  -Wl,-keep_dwarf_unwind -o ${BINDIR}/${EXECUTABLE} ${OBJECTS} ${LDFLAGS} ${PETSC_VEC_LIB} ${SLEPC_LIB}
