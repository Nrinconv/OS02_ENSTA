include Make_linux.inc
# include Make_msys2.inc
# include Make_osx.inc

# Utiliser le compilateur MPI
CXX = mpic++

# Activer OpenMP avec -fopenmp
CXXFLAGS = -std=c++17 -Wall -fopenmp
ifdef DEBUG
	CXXFLAGS += -g -O0 -fbounds-check -pedantic -D_GLIBCXX_DEBUG -fsanitize=address
else
	CXXFLAGS += -O3 -march=native
endif

# Bibliothèques à l'édition de liens : MPI, OpenMP, pthread et SDL2
LIB += -lmpi -lpthread -lSDL2

ALL = simulation.exe

default: help

all: $(ALL)

clean:
	@rm -fr *.o *.exe *~ *.d

.cpp.o:
	$(CXX) $(CXXFLAGS) -MMD -c $< -o $@

simulation.exe: display.o model.o simulation.o
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LIB)

help:
	@echo "Available targets : "
	@echo "    all            : compile all executables"
	@echo "Add DEBUG=yes to compile in debug"
	@echo "Configuration :"
	@echo "    CXX      :    $(CXX)"
	@echo "    CXXFLAGS :    $(CXXFLAGS)"

%.html: %.md
	pandoc -s --toc $< --css=./github-pandoc.css --metadata pagetitle="OS202 - TD1" -o $@

-include *.d
