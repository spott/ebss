include ../makefile.include

# To generate .clang_complete:
# make make PCC='~/.vim/bundle/clang_complete/bin/cc_args.py /opt/local/petsc-3.3-gcc-debug/bin/mpicxx' -B <target>

SOURCES=nonlinear_index.cpp # mat_diff.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=nonlinear
DEFAULT=all

all: $(SOURCES) $(EXECUTABLE)

mat_diff: mat_diff.o chkopts
	-${CLINKER} -o ${BINDIR}/mat_diff mat_diff.o ${LDFLAGS} ${PETSC_VEC_LIB} ${SLEPC_LIB}

mat_view: mat_examine.o chkopts
	-${CLINKER} -o ${BINDIR}/mat_view mat_examine.o ${LDFLAGS} ${PETSC_VEC_LIB} ${SLEPC_LIB}

vec_view: vec_view.o chkopts
	-${CLINKER} -o ${BINDIR}/vec_view vec_view.o ${LDFLAGS} ${PETSC_VEC_LIB} ${SLEPC_LIB}

fourier: fourier.o chkopts
	#-${CLINKER} -o ${BINDIR}/fourier fourier.o ${LDFLAGS} ${PETSC_VEC_LIB} ${SLEPC_LIB}

#etok: etok.o chkopts
	#-${CLINKER} -o ${BINDIR}/etok etok.o ${LDFLAGS} ${PETSC_VEC_LIB} ${SLEPC_LIB}

#etok_computer: etok_computer.o chkopts
	#-${CLINKER} -o ${BINDIR}/etok_computer etok_computer.o ${LDFLAGS} ${PETSC_VEC_LIB} ${SLEPC_LIB}
#LDFLAGS+=-lboost_mpi-rt -lboost_serialization-rt

#nonlinear-index: nonlinear_index.o chkopts
	#-${CLINKER} -o ${BINDIR}/nonlinear_index nonlinear_index.o ${LDFLAGS} ${PETSC_VEC_LIB} ${SLEPC_LIB}

${EXECUTABLE}: ${OBJECTS}  chkopts
	-${CLINKER} -o ${BINDIR}/${EXECUTABLE} ${OBJECTS} ${LDFLAGS} ${PETSC_VEC_LIB} ${SLEPC_LIB}

#.phony clean: 
	#rm *.o
.PHONEY: clean
clean::
	rm *.o
