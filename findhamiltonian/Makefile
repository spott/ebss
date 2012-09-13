include ../makefile.include

SOURCES=findhamiltonian.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=findhamiltonian
DEFAULT=all

all: $(SOURCES) $(EXECUTABLE)

${EXECUTABLE}: ${OBJECTS}  chkopts
	-${CLINKER} -o ${BINDIR}/${EXECUTABLE} ${OBJECTS} ${LDFLAGS} ${PETSC_VEC_LIB} ${SLEPC_LIB}
