
all:
	$(MAKE) -C src/vpp
	$(MAKE) -C src/vpp-mex
	$(MAKE) -C src/1_preprocessing/ponomarenko
	$(MAKE) -C src/2_stabilization/estadeo
	$(MAKE) -C src/3_oflow/tvl1flow_3
	$(MAKE) -C src/4_denoising/kalman
	$(MAKE) -C src/5_deblurring/fba-ipol

clean:
	$(MAKE) clean -C src/vpp
	$(MAKE) clean -C src/vpp-mex
	$(MAKE) clean -C src/1_preprocessing/ponomarenko
	$(MAKE) clean -C src/2_stabilization/estadeo
	$(MAKE) clean -C src/3_oflow/tvl1flow_3
	$(MAKE) clean -C src/4_denoising/kalman
	$(MAKE) clean -C src/5_deblurring/fba-ipol


