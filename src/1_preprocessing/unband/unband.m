#!/usr/bin/env octave-cli

% get list of arguments passed to script
arg_list = argv();
inputName = arg_list{1};
outputName = arg_list{2};

% load required Octave packages
pkg load statistics
pkg load image

% add vpp to the path
addpath([fileparts(mfilename('fullpath')) '/../../vpp-mex']);


% parameters
rfactor = 6;
filt_sz = 9;
sigma_x = 5;

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

v = vpp_read_frame(input);
while ~isscalar(v)
	% helper matrices
	mv = ones(size(v,1)-20+1,1);
	mh = ones(1,size(v,2));

	tmp = v(20:end,:);
	%tmp = tmp - mean(tmp(:)) + mmm;

	% vertical mean and std of temporal median
	%vert_mean = median(tmp);
	%vert_std = 1.4826 * mad(tmp,1);
	vert_mean = mean(tmp);
	vert_std = std(tmp);

	% smoothing with gaussian filter
	gvert_mean = imfilter(vert_mean, fspecial('gaussian',[1 11], 3), 'symmetric');
	gvert_std  = imfilter(vert_std , fspecial('gaussian',[1 11], 3), 'symmetric');
		
	% smoothing with joint bilateral filter
	signal1 = vert_mean;
	signal2 = vert_std;
	out1 = zeros(size(signal1));
	out2 = zeros(size(signal2));
	sigma_r1 = rfactor * mad(vert_mean - gvert_mean,1)*1.4826;
	sigma_r2 = rfactor * mad(vert_std  - gvert_std,1) *1.4826;
	%disp([sigma_r1 sigma_r2])
	signal1_bc = [signal1(1)*ones(1,filt_sz) signal1 signal1(end)*ones(1,filt_sz)];
	signal2_bc = [signal2(1)*ones(1,filt_sz) signal2 signal2(end)*ones(1,filt_sz)];
	filt_x = fspecial('gaussian',[1 2*filt_sz+1], 3);
	for j = filt_sz + [1:length(signal1)],
		vals1 = signal1_bc(j + [-filt_sz:filt_sz]);
		vals2 = signal2_bc(j + [-filt_sz:filt_sz]);
		filt_r1 = exp(-1/2/sigma_r1^2*(vals1 - signal1_bc(j)).^2);
		filt_r2 = exp(-1/2/sigma_r2^2*(vals2 - signal2_bc(j)).^2);
		filt_w = filt_r1.*filt_r2.*filt_x;
		out1(j - filt_sz) = sum(filt_w.*vals1)/sum(filt_w);
		out2(j - filt_sz) = sum(filt_w.*vals2)/sum(filt_w);
	end
	svert_mean = out1;
	svert_std  = out2;


	%clf, plot(vert_mean, 'r' ), hold on, plot(svert_mean,'k'), hold off
	%clf, plot(vert_std, 'r' ), hold on, plot(svert_std,'k'), hold off
	%ylim([8880 8920])
	%pause(.1)

	tmp = (tmp - mv*vert_mean).*(mv*(svert_std./vert_std)) + mv*svert_mean;

	% horizontal mean and std of temporal median
	%horz_mean = median(tmp,2);
	%horz_std = mad(tmp,1,2);
	horz_mean = mean(tmp,2);
	horz_std = std(tmp,0,2);

	% smoothing with gaussian kernel
	ghorz_mean = imfilter(horz_mean, fspecial('gaussian',[10 1], 3), 'symmetric');
	ghorz_std  = imfilter(horz_std , fspecial('gaussian',[10 1], 3), 'symmetric');

	% smoothing with joint bilateral filter
	signal1 = horz_mean';
	signal2 = horz_std';
	out1 = zeros(size(signal1));
	out2 = zeros(size(signal2));
	sigma_r1 = rfactor * mad(horz_mean - ghorz_mean,1)*1.4826;
	sigma_r2 = rfactor * mad(horz_std  - ghorz_std,1) *1.4826;
	signal1_bc = [signal1(1)*ones(1,filt_sz) signal1 signal1(end)*ones(1,filt_sz)];
	signal2_bc = [signal2(1)*ones(1,filt_sz) signal2 signal2(end)*ones(1,filt_sz)];
	filt_x = fspecial('gaussian',[1 2*filt_sz+1], 3);
	for j = filt_sz + [1:length(signal1)],
		vals1 = signal1_bc(j + [-filt_sz:filt_sz]);
		vals2 = signal2_bc(j + [-filt_sz:filt_sz]);
		filt_r1 = exp(-1/2/sigma_r1^2*(vals1 - signal1_bc(j)).^2);
		filt_r2 = exp(-1/2/sigma_r2^2*(vals2 - signal2_bc(j)).^2);
		filt_w = filt_r1.*filt_r2.*filt_x;
		out1(j - filt_sz) = sum(filt_w.*vals1)/sum(filt_w);
		out2(j - filt_sz) = sum(filt_w.*vals2)/sum(filt_w);
	end
	shorz_mean = out1';
	shorz_std  = out2';

	%clf, plot(svert_mean - vert_mean, 'r' ), hold off
	%ylim([-2 2])
	%clf, plot(horz_mean, 'r' ), hold on, plot(shorz_mean,'k'), hold off
	%ylim([8880 8920])
	%pause(.1)

	uv = v;
	uv(20:end,:) = (tmp - horz_mean*mh).*((shorz_std./horz_std)*mh) + shorz_mean*mh;

	if ~vpp_write_frame(output, uv)
		break
	end
	v = vpp_read_frame(input);
end
vpp_close(input);
vpp_close(output);

