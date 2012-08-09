include ../makefile.include

SOURCES=test_parameters.cpp ${TDSEBASE}/common/parameters.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=test
DEFAULT=all

all: $(SOURCES) $(EXECUTABLE)

${EXECUTABLE}: ${OBJECTS}  chkopts 
	-${CLINKER} -o ${BINDIR}/${EXECUTABLE} ${OBJECTS} ${LDFLAGS} ${PETSC_VEC_LIB}

