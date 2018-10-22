CC=c99
OMPFLAGS=-fopenmp
CFLAGS=-O3 -I../../vpp/ $(OMPFLAGS)

VPP_OBJ=../../vpp/vpp.o
OBJ=main.o tvl1flow_lib.o bicubic_interpolation.o mask.o zoom.o xmalloc.o

tvl1flow: $(VPP_OBJ) $(OBJ) 
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS) -ltiff -ljpeg -lpng -lm $(OMPFLAGS)


clean:
	rm -f $(OBJ) tvl1flow
