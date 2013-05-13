include ../../../makefile.include

SOURCES=examples.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=examples
DEFAULT=all

all: $(SOURCES) $(EXECUTABLE)

${EXECUTABLE}: ${OBJECTS}  chkopts 
	-${CLINKER} -o ${BINDIR}/${EXECUTABLE} ${OBJECTS} ${LDFLAGS} ${PETSC_VEC_LIB}

