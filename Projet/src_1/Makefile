include Make_linux.inc
#include Make_msys2.inc
#include Make_osx.inc

CXXFLAGS = -std=c++17 -Wall -O3 -march=native -fopenmp

ALL = simulation_parallel.exe 

default: help

all: $(ALL)

clean:
	@rm -fr *.o *.exe *~ *.d

.cpp.o:
	$(CXX) $(CXXFLAGS) -MMD -c $< -o $@

simulation_parallel.exe: display.o model.o simulation.o
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LIB)	

help:
	@echo "Available targets : "
	@echo "    all            : compile parallel version"
	@echo "Configuration :"
	@echo "    CXX      :    $(CXX)"
	@echo "    CXXFLAGS :    $(CXXFLAGS)"

%.html: %.md
	pandoc -s --toc $< --css=./github-pandoc.css --metadata pagetitle="OS202 - TD1" -o $@

-include *.d
