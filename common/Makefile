include ../makefile.include

SOURCES=parameters_test.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=test
DEFAULT=all

all: $(SOURCES) $(EXECUTABLE)

${EXECUTABLE}: ${OBJECTS}  chkopts 
	-${CLINKER} -o ${BINDIR}/${EXECUTABLE} ${OBJECTS} ${LDFLAGS} ${PETSC_VEC_LIB}

