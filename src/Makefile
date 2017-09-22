#
# Definition of MACROS

BINROOT=./
CXX = g++
CXXFLAGS = -Wall
LIBS = -lm
AR = ar rcs
LINKFLAGS = $(CXXFLAGS)

EX = ChemNetworks-2.2.exe
EXLIB = libChemNetworks-2.2.a

ARCH?= gcc
include arch/$(ARCH)

SRC = $(wildcard *.cpp)
INC = $(wildcard *.h)
OBJ = $(filter-out main.o, $(SRC:.cpp=.o))

#
# -- target : 	Dependencies
# --		Rule to create target

$(EX): main.o $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(LINKFLAGS) main.o $(OBJ) $(LIBS)

####################################################################

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

lib:
	$(AR) $(EXLIB) $(OBJ)

####################################################################

clean: 
	rm -rf *.o $(EX) $(EX).dSYM $(EXLIB)

#
# -- Simple dependencies
main.o : $(INC)
$(OBJ) : $(INC)
