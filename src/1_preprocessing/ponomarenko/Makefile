# Makefile created with Create Makefile.sh 13/12/2011

COPT      = -O99 -g
CFLAGS    = $(COPT) -fopenmp
CSTRICT   = -Wall -Wextra -ansi
CXXOPT    = -O99 -g
CXXFLAGS  = $(CXXOPT) -fopenmp
CXXSTRICT = -Wall -Wextra -ansi
LDFLAGS   = -lpng -ljpeg -ltiff -lm -lfftw3f -lfftw3f_threads -lgomp
EXEC      = ponomarenko

OBJ=./iio/iio.o\
	 ./framework/CFramework.o\
	 ./framework/libparser.o\
	 ./framework/CHistogram.o\
	 ./framework/CImage.o\
	 ./framework/CInspector.o\
	 ./framework/operations.o\
	 ./curve_filter.o\
	 ./algo.o\
	 ./main.o


default: $(EXEC)
all: $(EXEC)

# ------- C files -------
# ------- C++ files -------
./framework/CFramework.o: ./framework/CFramework.cpp ./framework/CFramework.h
	$(CXX) $(CXXFLAGS) -c ./framework/CFramework.cpp -o ./framework/CFramework.o

./framework/libparser.o: ./framework/libparser.cpp ./framework/libparser.h
	$(CXX) $(CXXFLAGS) -c ./framework/libparser.cpp -o ./framework/libparser.o

./framework/CHistogram.o: ./framework/CHistogram.cpp
	$(CXX) $(CXXFLAGS) -c ./framework/CHistogram.cpp -o ./framework/CHistogram.o

./framework/CImage.o: ./framework/CImage.cpp ./framework/CImage.h
	$(CXX) $(CXXFLAGS) -Iiio -c ./framework/CImage.cpp -o ./framework/CImage.o

./framework/CInspector.o: ./framework/CInspector.cpp ./framework/CInspector.h
	$(CXX) $(CXXFLAGS) -c ./framework/CInspector.cpp -o ./framework/CInspector.o

./framework/operations.o: ./framework/operations.cpp
	$(CXX) $(CXXFLAGS) -c ./framework/operations.cpp -o ./framework/operations.o

./curve_filter.o: ./curve_filter.cpp ./curve_filter.h
	$(CXX) $(CXXFLAGS) -c ./curve_filter.cpp -o ./curve_filter.o

./algo.o: ./algo.cpp ./algo.h
	$(CXX) $(CXXFLAGS) -c ./algo.cpp -o ./algo.o

./main.o: ./main.cpp
	$(CXX) $(CXXFLAGS) -c ./main.cpp -o ./main.o


# ------- Main -------
$(EXEC): $(OBJ)
	$(CXX) $(OBJ) $(LDFLAGS) -o $(EXEC)

lint:
	$(MAKE) CFLAGS="$(CFLAGS) $(CSTRICT)" CXXFLAGS="$(CXXFLAGS) $(CXXSTRICT)"

clean:
	rm -f $(OBJ)

distclean: clean
	rm -f $(EXEC)

