CFLAGS ?= -march=native -O3 -Wall -Wextra
CFLAGS += -I. -I../../vpp/ -Irepo/src/ -Irepo/lib -g -std=gnu99
LDLIBS = -lm -fopenmp

OBJ = ../../vpp/vpp.o
BIN = kalman

default: $(BIN)

kalman: main.o $(OBJ) repo/lib/argparse/argparse.o
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS)  -ltiff -ljpeg -lpng -lm -lfftw3f

clean: ; @$(RM) $(OBJ) *.o repo/lib/argparse/*.o
dist_clean: ; @$(RM) $(BIN) $(OBJ) *.o repo/lib/argparse/*.o
.PRECIOUS: %.o

