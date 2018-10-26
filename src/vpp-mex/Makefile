
BINS=vpp_init_input.mex\
	  vpp_init_output.mex\
	  vpp_read_frame.mex\
	  vpp_write_frame.mex\
	  vpp_close.mex

default: $(BINS)

%.mex: %.c
	mkoctfile --mex $^ vpp/vpp.c -Ivpp

clean:
	-rm $(BINS)
	-rm *.o
