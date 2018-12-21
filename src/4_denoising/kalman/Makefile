CFLAGS ?= -march=native -O3 -Wall -Wextra
CFLAGS += -I. -I../../vpp/ -Irepo/src/ -Irepo/lib -g -std=gnu99
LDLIBS = -lm -fopenmp

BIN = kalman

default: $(BIN)

kalman: main.o repo/lib/argparse/argparse.o
	$(CC) $(LDFLAGS) -o $@ $^ ../../vpp/vpp.o $(LDLIBS)  -ltiff -ljpeg -lpng -lm -lfftw3f

clean: ; @$(RM) *.o repo/lib/argparse/*.o
dist_clean: ; @$(RM) $(BIN) *.o repo/lib/argparse/*.o
.PRECIOUS: %.o

