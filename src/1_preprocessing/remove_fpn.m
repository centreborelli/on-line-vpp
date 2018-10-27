#!/usr/bin/env octave-cli

% get list of arguments passed to script
arg_list = argv();
inputName = arg_list{1};
outputName = arg_list{2};

% load required Octave packages
pkg load statistics
pkg load image

% add vpp to the path
addpath([fileparts(mfilename('fullpath')) '/../vpp-mex']);

% parameters
blursigma = 5;
lowpass = fspecial('gaussian', 4*blursigma+1, blursigma);

input = vpp_init_input(inputName);
if ~input
	disp('cannot load input');
	exit(1);
end

output = vpp_init_output(outputName, input);
if ~output
	disp('cannot load output');
	exit(1);
end

frames = zeros(input(1), input(2), input(3), 0);

v = vpp_read_frame(input);
while v
	frames(:,:,:,end+1) = v;
	v = vpp_read_frame(input);
end

med = median(frames, 4);
fpn = med - imfilter(med, lowpass, 'same', 'symmetric');
frames = frames - fpn;

for f=1:size(frames,4)
	if ~vpp_write_frame(output, frames(:,:,:,f))
		break
	end
end
vpp_close(input);
vpp_close(output);

