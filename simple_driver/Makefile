#
# Definition of MACROS

BINROOT=./
CXX = g++
CXXFLAGS = -Wall
LIBS = -lm
AR = ar rcs

EX = simple_driver

ARCH?= gcc
include arch/$(ARCH)

SRC = $(wildcard *.cpp)
INC = $(wildcard *.h)
OBJ = $(SRC:.cpp=.o)

#
# -- target : 	Dependencies
# --		Rule to create target

$(EX): $(OBJ)
	$(CXX) -o $@ $(CXXFLAGS) $(OBJ) $(LIBS)

####################################################################

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

lib:
	$(AR) $(EXLIB) $(OBJ)

####################################################################

clean: 
	rm -f *.o $(EX)

#
# -- Simple dependencies
$(OBJ) : $(INC)
