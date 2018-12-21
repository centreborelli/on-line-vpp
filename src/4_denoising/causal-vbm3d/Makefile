# C++ source code
CXXSRC	= main.cpp \
		vbm3d.cpp \
		../../vpp/vpp.c \
		Utilities.cpp \
		lib_transforms.cpp

# all source code
SRC	= $(CSRC) $(CXXSRC)

# C objects
COBJ	= $(CSRC:.c=.o)
# C++ objects
CXXOBJ	= $(CXXSRC:.cpp=.o)
# all objects
OBJ	= $(COBJ) $(CXXOBJ)
# binary target
BIN	= vbm3d

# C optimization flags
COPT	= -O3 -ftree-vectorize -funroll-loops
#COPT	= 

# C++ optimization flags
CXXOPT	= $(COPT)

# C compilation flags
CFLAGS	= $(COPT) -I../../vpp -Wno-write-strings -ansi -std=c99# -g
# C++ compilation flags
CXXFLAGS	= $(CXXOPT)	-I../../vpp -Wno-write-strings -Wno-deprecated -ansi -std=c++11# -g -fsanitize=address
# link flags
LDFLAGS	= -lpng -lm -lfftw3f -ltiff -ljpeg# -fsanitize=address

# use openMP with `make OMP=1`
ifdef OMP
CFLAGS	+= -fopenmp
CXXFLAGS	+= -fopenmp
LDFLAGS += -lgomp
else
CFLAGS	+= -Wno-unknown-pragmas
CXXFLAGS  += -Wno-unknown-pragmas
endif

# partial compilation of C source code
%.o: %.c %.h
	$(CC) -c -o $@  $< $(CFLAGS)
# partial compilation of C++ source code
%.o: %.cpp %.h
	$(CXX) -c -o $@  $< $(CXXFLAGS)

# link all the object code
$(BIN): $(OBJ) $(LIBDEPS)
	$(CXX) -o $@ $(OBJ) $(LDFLAGS)

clean:
	rm *.o $(BIN)
