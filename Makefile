
all: vpp tvl1flow kalman

vpp:
	$(MAKE) -C src/$@

tvl1flow:
	$(MAKE) -C src/3_oflow/tvl1flow_3

kalman:
	$(MAKE) -C src/4_denoising/kalman
