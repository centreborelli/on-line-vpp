#! /bin/bash
# Script implementing the following pipeline
#    add noise and create duplicate
#    compute flow and create duplicate
#    detect occlusions
#    denoise
#    write output
#
# readvid ─🠆 add noise ─🠆 dup ─🠆 tvl1flow ─🠆 dup ─🠆 occlusions ─🠆 kalman ─🠆 writevid
#                          │                  │                    🡑 🡑
#                          │                  └────(oflow)─────────┘ │
#                          └─────(noisy)─────────────────────────────┘

sigma=$2;


# read input
mkfifo noisy
mkfifo oflow
bin/readvid "$1" - \
	| bin/vlambda - "randn $sigma * x +" \
	| bin/vp dup - noisy - \
	| bin/tvl1flow - - \
	| bin/vp dup - oflow - \
	| bin/vlambda - "x(0,0)[0] x(-1,0)[0] - x(0,0)[1] x(0,-1)[1] - + 0.5 > 255 *" -o - \
	| bin/kalman -i noisy -o oflow -k - -s $sigma -d - \
#	| bin/fba 5 128 11 3 - - \
	| bin/writevid - /tmp/d%03d.png
rm noisy oflow

