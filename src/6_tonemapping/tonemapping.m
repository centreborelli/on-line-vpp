#!/usr/bin/env octave-cli

%%% Warning:
% Messages must be written in stderr, because stdout is used for video flux.

% Get list of arguments passed to script
arg_list = argv();
if nargin < 2 || nargin > 3
    error(['Incorrect number of arguments.\n' ...
        'Usage: tonemapping.m input output [ratioTarget]\n']);
end
inputName = arg_list{1}; % can be '-' for stdin
outputName = arg_list{2}; % can be '-' for stdout
if nargin == 3
    ratioTarget = arg_list{3};
else
    ratioTarget = 1/4;
end

%%% Add vpp mex files to path
addpath([fileparts(mfilename('fullpath')) '/../vpp-mex']);
addpath(fileparts(mfilename('fullpath')));

%%% Load required Octave packages
pkg load statistics
pkg load image

%%% Input
inputHandle = vpp_init_input(inputName);
if ~inputHandle, error('cannot load input'); end

%%% Define size of output
H = inputHandle(1); % height of input
W = inputHandle(2); % width of input
D = 3;              % depth: output has color

%%% Output
outputHandle = vpp_init_output(outputName, [H,W,D]);
if ~outputHandle, error('cannot load output'); end

%%% Get colormap (create variable "cmap")
warning ('off', 'Octave:data-file-in-path');
load([fileparts(mfilename('fullpath')) '/colormap.mat']);
if ~exist('cmap','var'), error('Couldn''t load the colormap.'); end
cmap = single(cmap);

%%% MGF parameters
epsSqrt = 20;   % initialize epsilon. Its value is set automatically afterwards
radius  = 1;    % small images, need to avoid contrast halo
nIt     = 1;    % for speed
nScales = 0;    % filter all scales
gamma   = 1;    % evolution of gamma over the scales: eps <-- eps / scale^gamma
                %   - with gamma>0, reduce eps for low scales
                %   - with gamma=0, keep eps constant

%%% Estimation of epsilon: parameters and initializations
ratioThreshold = .05;
step = NaN;

%%% Get first frame; then loop until there is no frame available
u = vpp_read_frame(inputHandle); % single
n = 1;
while ~isscalar(u)

    %%% Place input in a range where squared values does not explode
    %%% MGF does not need the input to be in a specific range, yet large values
    %%% cause numerical errors in single precision.
    %%% Use u = (u-mean(u(:)))/std(u(:)) if needed, or use double.
    u = u - mean(u(:));

    %%% compute only statistics on valid pixels
    %%% u=0 happens due to the stabilization
    mask = u~=0;

    %%% Estimation of epsilon: initializations
    uStd  = mad(u(mask),0); % std(u(:));
    if uStd==0, uStd = eps; end % avoid ratio = NaN
    ratio = inf;
    nIter = 1;

    %%% Compute parameter epsilon, and base + detail decomposition
    while abs(ratio - ratioTarget) > ratioThreshold

        %%% Compute decomposition
        b = MGFoctave(u, u, epsSqrt^2, radius, nIt, nScales, gamma);
        d = u - b;

        %%% Compute ratio std(detail)/std(image)
        dStd  = mad(d(mask),0); % std(d(:));
        ratio = dStd / uStd;

        fprintf(2,...
            'Img #%04d\tIter #%d:\tEpsilon = %f\tRatio = %f\tStep = %f\n',...
            n, nIter, epsSqrt, ratio, step);

        %%% Update values
        epsSqrtPrev = epsSqrt;
        epsSqrt     = epsSqrt * ratioTarget / ratio;
        if nIter == 1 % handle brutal changes
            epsSqrt = .25*epsSqrt + .75*epsSqrtPrev; % temporal smoothing
        end
        step        = epsSqrt - epsSqrtPrev;
        nIter       = nIter + 1;

        %%% Saturate epsSqrt
        %%% For mostly flat images (no detail), ratioTarget can't be reached
        if epsSqrt > 200
            epsSqrt = 200;
            if epsSqrtPrev == 200, break; end
        end

        %%% If estimation does not converge
        %%% This should not happen because of the condition above
        if nIter > 15, break; end
    end

    %%% Compress the base layer (in [-.5,+.5])
    bMed = median(b(:));
    bMad = median(abs(u(:) - bMed)); % same as mad(b(:),1);
    b    = max(-.5,min(+.5, (b-bMed)./(36*bMad) ));

    %%% Compress the detail layer (in [-.5,+.5])
    % dMad = mad(d(:),1);
    % d    = d ./ (12*dMad);
    if n == 1
        dMin = min(d(:));
        dMax = max(d(:));
    else
        dMin = .25*min(d(:)) + .75*dMin;
        dMax = .25*max(d(:)) + .75*dMax;
    end
    d = (d-dMin)./(dMax-dMin) - .5;

    %%% Gray tone-mapped image
    %%% We prefer clipping in black than in white, hence the .4 instead of .5
    v = .4 + b/2 + d;

    %%% Count black and white clipping
    vClipBlack = sum( v(:) < 0 ) / numel(v);
    vClipWhite = sum( v(:) > 1 ) / numel(v);
    fprintf(2,...
        'Img #%04d\tClipping: %.2f%% in black and %.2f%% in white.\n',...
        n, vClipBlack*100, vClipWhite*100);

    %%% Clip
    v = max(0,min(1, v ));

    %%% Compute color coefficients using the (normalized+clipped) input image
    u = max(0,min(1, .5 + (u-bMed)./(12*bMad) ));
    u = single(ind2rgb(uint16( u*(2^16-1) ), cmap));
    c = u ./ (sum(u,3)/3);

    %%% Apply color coefficient to vGray to get final colored tonemapped image
    v = 255 * v .* c;

    %%% Write frame
    if ~vpp_write_frame(outputHandle, v), break; end

    %%% Read next frame
    u = vpp_read_frame(inputHandle);
    n = n+1;
end

vpp_close(inputHandle);
vpp_close(outputHandle);
