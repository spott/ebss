include ../makefile.include

SOURCES=propagate.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=propagate
DEFAULT=all

all: $(SOURCES) $(EXECUTABLE)

${EXECUTABLE}: ${OBJECTS}  chkopts
	-${CLINKER} -o ${BINDIR}/${EXECUTABLE} ${OBJECTS} ${LDFLAGS} ${PETSC_VEC_LIB} ${SLEPC_LIB}
