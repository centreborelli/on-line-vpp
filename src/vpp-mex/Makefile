BINS=vpp_init_input.mex\
	  vpp_init_output.mex\
	  vpp_read_frame.mex\
	  vpp_write_frame.mex\
	  vpp_close.mex

default: $(BINS)

%.mex: %.c
	mkoctfile ${CFLAGS} --mex $^ ../vpp/vpp.o -I../vpp

clean:
	rm -f $(BINS)

