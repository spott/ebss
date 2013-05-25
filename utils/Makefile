include ../makefile.include

# To generate .clang_complete:
# make make PCC='~/.vim/bundle/clang_complete/bin/cc_args.py /opt/local/petsc-3.3-gcc-debug/bin/mpicxx' -B <target>

SOURCES=etok_computer.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=etok_computer
DEFAULT=all

all: $(SOURCES) $(EXECUTABLE)

fourier: fourier.o chkopts
	-${CLINKER} -o ${BINDIR}/fourier fourier.o ${LDFLAGS} ${PETSC_VEC_LIB} ${SLEPC_LIB}

etok: etok.o chkopts
	-${CLINKER} -o ${BINDIR}/etok etok.o ${LDFLAGS} ${PETSC_VEC_LIB} ${SLEPC_LIB}

etok_computer: etok_computer.o chkopts
	-${CLINKER} -o ${BINDIR}/etok_computer etok_computer.o ${LDFLAGS} ${PETSC_VEC_LIB} ${SLEPC_LIB}

#${EXECUTABLE}: ${OBJECTS}  chkopts
#	-${CLINKER} -o ${BINDIR}/${EXECUTABLE} ${OBJECTS} ${LDFLAGS} ${PETSC_VEC_LIB} ${SLEPC_LIB}

#.phony clean: 
	#rm *.o
.PHONEY: clean
clean::
	rm *.o
