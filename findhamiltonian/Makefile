include ../makefile.include

TESTS=bspline/bspline-test.cpp
SOURCES=findbasis.cpp
OBJECTS=$(SOURCES:.cpp=.o)
TOBJECTS=$(TESTS:.cpp=.o)
EXECUTABLE=findbasis
TEXECUTABLE=bspline/bspline-test
DEFAULT=all

all: $(SOURCES) $(EXECUTABLE)

${TEXECUTABLE}: ${TOBJECTS}  chkopts
	-${CLINKER} -o ${TEXECUTABLE} ${TOBJECTS} ${LDFLAGS} ${PETSC_VEC_LIB} ${SLEPC_LIB}

${EXECUTABLE}: ${OBJECTS}  chkopts
	-${CLINKER} -o ${BINDIR}/${EXECUTABLE} ${OBJECTS} ${LDFLAGS} ${PETSC_VEC_LIB} ${SLEPC_LIB}
