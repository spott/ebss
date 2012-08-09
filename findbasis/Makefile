include ../makefile.include

TESTS=bspline/bspline-test.cpp ${TDSEBASE}/common/parameters.cpp
SOURCES=findbasis.cpp ${TDSEBASE}/common/parameters.cpp
OBJECTS=$(SOURCES:.cpp=.o)
TOBJECTS=$(TESTS:.cpp=.o)
EXECUTABLE=findbasis
TEXECUTABLE=bspline/bspline-test
DEFAULT=all

all: $(SOURCES) $(EXECUTABLE)

${TEXECUTABLE}: ${TOBJECTS}  chkopts
	-${CLINKER} -o ${TEXECUTABLE} ${TOBJECTS} ${LDFLAGS} ${PETSC_VEC_LIB}

${EXECUTABLE}: ${OBJECTS}  chkopts
	-${CLINKER} -o ${BINDIR}/${EXECUTABLE} ${OBJECTS} ${LDFLAGS} ${PETSC_VEC_LIB}
