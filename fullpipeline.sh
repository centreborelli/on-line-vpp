#!/usr/bin/env bash

noisy=noisy
oflow=oflow
pono=pono
sigmas=sigmas

rm $noisy $oflow $pono $sigmas
mkfifo $noisy
mkfifo $oflow
mkfifo $pono
mkfifo $sigmas

# optical flow parameters
DW=0.1 # weight of data attachment term 0.1 ~ very smooth 0.2 ~ noisy
FS=1   # finest scale (0 image scale, 1 one coarser level, 2 more coarse, etc...
oflow_params="0 0.25 $DW 0.3 100 $FS 0.5 5 0.01 0" 

bin/readvid "$1" - \
	`# preprocessing` \
	| ./src/1_preprocessing/remove_fpn.m - - \
	| ./src/1_preprocessing/unband/unband.m - - \
	| ./bin/vp dup - - $pono \
	`# stabilization` \
	| bin/estadeo - -o - \
	`# nlkalman (with opticalflow and occlusion detection)` \
	| bin/vp dup - $noisy - \
	| bin/tvl1flow - - $oflow_params \
	| bin/vp dup - $oflow - \
	| bin/vlambda - "x(0,0)[0] x(-1,0)[0] - x(0,0)[1] x(0,-1)[1] - + 0.5 > 255 *" -o - \
	| bin/kalman -i $noisy -o $oflow -k - -s $sigmas -d - -v \
	`# tonemapping` \
	| ./src/6_tonemapping/tonemapping.m - - \
	`# save` \
	| bin/writevid - $2 \
	`# estimate noise` \
	| bin/ponomarenko -b 1 $pono $sigmas

rm $noisy $oflow $pono $sigmas

