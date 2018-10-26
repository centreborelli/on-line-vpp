
all:
	$(MAKE) -C src/vpp
	$(MAKE) -C src/vpp-mex
	$(MAKE) -C src/1_preprocessing/ponomarenko
	$(MAKE) -C src/3_oflow/tvl1flow_3
	$(MAKE) -C src/4_denoising/kalman
	$(MAKE) -C src/5_deblurring/fba-ipol

