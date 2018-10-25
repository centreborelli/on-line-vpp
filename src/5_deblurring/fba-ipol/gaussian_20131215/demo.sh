#! /bin/sh
# Demo script for gaussian_demo
# Pascal Getreuer 2013

# Echo shell commands
set -v

# Set parameters
sigma=5
algo=deriche
K=3
tol=1e-6

# Perform Gaussian convolution on the image einstein.png
./gaussian_demo -s $sigma -a $algo -K $K -t $tol einstein.png blurred.png
