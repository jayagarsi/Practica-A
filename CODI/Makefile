TARGET = basics greedy local_search metaheuristic
CXXFLAGS = -ansi -O3 -fpermissive -std=c++17 
OBJS = ./HEADERS/Random.o ./HEADERS/Timer.o
CPLOBJS = ./HEADERS/Random.o ./HEADERS/Timer.o

SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic
GCC = gcc
CCC = g++
CCOPT = -m64 -O -fPIC -fexceptions -DNDEBUG -DIL_STD -std=c++17 -fpermissive -w

CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR)

all: ${TARGET}

basics: basics.cpp ${OBJ}
	${CCC} ${CXXFLAGS} -o $@ $^

greedy: greedy.cpp $(OBJS)
	${CCC} ${CXXFLAGS} -o $@ $^

local_search: local_search.cpp $(OBJS)
	${CCC} ${CXXFLAGS} -o $@ $^

metaheuristic: metaheuristic.cpp $(OBJS)
	${CCC} ${CXXFLAGS} -o $@ $^

clean:
	@rm -f *~ *.o ${TARGET} core


