include ../makefile.include

# To generate .clang_complete:
# make make PCC='~/.vim/bundle/clang_complete/bin/cc_args.py /opt/local/petsc-3.3-gcc-debug/bin/mpicxx' -B <target>

SOURCES=etok.cpp field_free_propagate.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=etok
DEFAULT=all

all: $(SOURCES) $(EXECUTABLE)

${EXECUTABLE}: ${OBJECTS}  chkopts
	-${CLINKER} -o ${BINDIR}/${EXECUTABLE} ${OBJECTS} ${LDFLAGS} ${PETSC_VEC_LIB} ${SLEPC_LIB}

field_free_propagate: field_free_propagate.o  chkopts
	-${CLINKER} -o ${BINDIR}/field_free_propagate ${OBJECTS} ${LDFLAGS} ${PETSC_VEC_LIB} ${SLEPC_LIB}
