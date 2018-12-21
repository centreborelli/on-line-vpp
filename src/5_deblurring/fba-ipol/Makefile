CFLAGS ?= -march=native -O3 -Wall -Wextra

OBJ=main.o fba.o consistent_registration.o utils.o
OBJ+=gaussian_20131215/basic.o gaussian_20131215/filter_util.o gaussian_20131215/invert_matrix.o gaussian_20131215/gaussian_short_conv.o gaussian_20131215/gaussian_conv_vyv.o
OBJ+=tvl1flow_3/tvl1flow_lib.o

CFLAGS+=-Igaussian_20131215/ -I../../vpp/
LDLIBS+=-lm -lfftw3

default: fba

fba: ${OBJ}
	${CC} $^ ../../vpp/vpp.o -o $@ ${CFLAGS} ${LDLIBS}

clean: ; @$(RM) $(OBJ)
dist_clean: ; @$(RM) $(BIN) $(OBJ)
.PRECIOUS: %.o

