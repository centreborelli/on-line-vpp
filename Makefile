
all: vpp tvl1flow kalman fba

vpp:
	$(MAKE) -C src/$@

tvl1flow:
	$(MAKE) -C src/3_oflow/tvl1flow_3

kalman:
	$(MAKE) -C src/4_denoising/kalman

fba:
	$(MAKE) -C src/5_deblurring/fba-ipol

