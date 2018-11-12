#!/usr/bin/env bash

src=$1
dst=$2
dsttmp=$(dirname $dst)_intermediate

OUTPUT_INTERMEDIATE=${OUTPUT_INTERMEDIATE:-1}
REMOVE_FPN=${REMOVE_FPN:-0}
STABILIZE=${STABILIZE:-1}
DENOISE=${DENOISE:-0}
DEBLUR=${DEBLUR:-0}

echo "Running full pipeline for sequence $src"
echo "		- fixed noise removal is $REMOVE_FPN"
echo "		- stabilization is $STABILIZE"
echo "		- debluring is $DEBLUR"

# create a temporary directory for the named pipes
tmp=$(mktemp -d)
trap "rm -rf $tmp" EXIT

# create named pipes used for the processing
function getfifo() {
	name=$1
	f=$tmp/$name
	mkfifo $f
	echo $f
}

noisy=$(getfifo noisy)
oflow=$(getfifo oflow)
pono=$(getfifo pono)
sigmas=$(getfifo sigmas)

# setup output files
if [ "$OUTPUT_INTERMEDIATE" -eq 1 ]; then
	mkdir -p $dsttmp 2>/dev/null
	outinput=$dsttmp/input.vpp
	outfpn=$dsttmp/fpn.vpp
	outunband=$dsttmp/unband.vpp
	outstab=$dsttmp/stab.vpp
	outdenoise=$dsttmp/denoise.vpp
	outdeblur=$dsttmp/deblur.vpp
	outtonemap=$dsttmp/tonemap.vpp
else
	outinput=/dev/null
	outfpn=/dev/null
	outunband=/dev/null
	outstab=/dev/null
	outdenoise=/dev/null
	outdeblur=/dev/null
	outtonemap=/dev/null
fi

# optical flow parameters
DW=0.1 # weight of data attachment term 0.1 ~ very smooth 0.2 ~ noisy
FS=1   # finest scale (0 image scale, 1 one coarser level, 2 more coarse, etc...
oflow_params="0 0.25 $DW 0.3 100 $FS 0.5 5 0.01 0" 

# define all the functions needed for the process
function remove_fpn() {
	./src/1_preprocessing/remove_fpn.m $1 - \
	| tee $outfpn
}

function unband() {
	./src/1_preprocessing/unband/unband.m $1 - \
	| tee $outunband
}

function stabilize() {
	bin/estadeo $1 -o - \
	| tee $outstab
}

function denoise() {
	bin/vp dup $1 $noisy - \
	| bin/tvl1flow - - $oflow_params \
	| bin/vp dup - $oflow - \
	| bin/vlambda - "x(0,0)[0] x(-1,0)[0] - x(0,0)[1] x(0,-1)[1] - + 0.5 > 255 *" -o - \
	| bin/kalman -i $noisy -o $oflow -k - -s $sigmas -d - -v \
	| tee $outdenoise
}

function deblur() {
	bin/fba 5 128 11 3 $1 - \
	| tee $outdeblur
}

function tonemap() {
	./src/6_tonemapping/tonemapping.m $1 - \
	| tee $outtonemap
}

# remove some of the features if needed
if [ "$REMOVE_FPN" -eq 0 ]; then
	function remove_fpn() {
		tee $outfpn
	}
fi
if [ "$STABILIZE" -eq 0 ]; then
	function stabilize() {
		tee $outstab
	}
fi
if [ "$DENOISE" -eq 0 ]; then
	function denoise() {
		tee $outdenoise
	}
fi
if [ "$DEBLUR" -eq 0 ]; then
	function deblur() {
		tee $outdeblur
	}
fi

# run everything!
bin/readvid "$src"'/*' - \
	| ./bin/vp dup - - $outinput \
	| remove_fpn - - \
	| unband - - \
	| ./bin/vp dup - - $pono \
	| stabilize - - \
	| denoise - - \
	| deblur - - \
	| tonemap - - \
	| bin/writevid - $dst \
	| bin/ponomarenko -b 1 $pono $sigmas \

