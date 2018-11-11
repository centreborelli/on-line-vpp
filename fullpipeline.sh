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

bin/readvid "$1" - \
	`# preprocessing` \
	| ./src/1_preprocessing/remove_fpn.m - - \
	| ./src/1_preprocessing/unband/unband.m - - \
	| ./bin/vp dup - - $pono \
	`# nlkalman (with opticalflow and occlusion detection)` \
	| bin/vp dup - $noisy - \
	| bin/tvl1flow - - \
	| bin/vp dup - $oflow - \
	| bin/vlambda - "x(0,0)[0] x(-1,0)[0] - x(0,0)[1] x(0,-1)[1] - + 0.5 > 255 *" -o - \
	| bin/kalman -i $noisy -o $oflow -k - -s $sigmas -d - -v \
	`# stabilization` \
	| bin/estadeo - -o - \
	`# tonemapping` \
	| ./src/6_tonemapping/tonemapping.m - - \
	`# save` \
	| bin/writevid - $2 \
	`# estimate noise` \
	| bin/ponomarenko -b 1 $pono $sigmas

rm $noisy $oflow $pono $sigmas

