include ../makefile.include

# To generate .clang_complete:
# make make PCC='~/.vim/bundle/clang_complete/bin/cc_args.py /opt/local/petsc-3.3-gcc-debug/bin/mpicxx' -B <target>

SOURCES=propagate.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=propagate
DEFAULT=all

all: $(SOURCES) $(EXECUTABLE)

${EXECUTABLE}: ${OBJECTS}  chkopts style chk_git
	-${CLINKER} -o ${BINDIR}/${EXECUTABLE} ${OBJECTS} ${LDFLAGS} ${PETSC_VEC_LIB} ${SLEPC_LIB}
