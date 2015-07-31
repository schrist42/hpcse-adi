 #	 
 # 	@author Fabian Thuering, 2015
 #  adapted
 #
 #	[DESCRIPTION]
 #	Makefile
 
CXX         = g++ #clang++

# === Sources ===
EXE         = main
SRC_CPP     = $(wildcard *.cpp)
OBJECTS     = $(patsubst %.cpp,objects/%.o,$(SRC_CPP))

# === Compiler Flags ===
WARNINGS    = -Wall
DEFINES     =  
DEBUG       =
INCLUDE     = 
OPT         = -O3 -march=native -g
CXXSTD      = -std=c++11
CXXFLAGS    = $(DEBUG) $(INCLUDE) $(CXXSTD) $(OPT) $(WARNINGS) $(DEFINES)
LDFLAGS     = -lboost_program_options -lboost_filesystem -lboost_system -lGL -lGLU -lglut

# === Compilation ===
.PHONY: all
all: $(EXE)

$(OBJECTS) : | objects

objects :
	mkdir objects


# Compile c++
objects/%.o : %.cpp
	$(CXX) -c $(CXXFLAGS) -o $@ $<
	
$(EXE) : $(OBJECTS)
	$(CXX) -o $(EXE) $(CXXFLAGS) $(OBJECTS) $(LDFLAGS)

# === Cleanup ===
.PHONY: clean
clean:	
	rm -f $(EXE) objects/*.o
