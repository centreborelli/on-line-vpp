#!/usr/bin/env octave-cli

%%% Warning:
% Messages must be written in stderr, because stdout is used for video flux.

% Get list of arguments passed to script
arg_list = argv();
if nargin < 2 || nargin > 4
    error(['Incorrect number of arguments.\n' ...
        'Usage: tonemapping.m input output [method [method-param]]\n']);
end
inputName = arg_list{1}; % can be '-' for stdin
outputName = arg_list{2}; % can be '-' for stdout
if nargin < 3
    method = 'localStd';
    increaseEps = 2;
else
    method = arg_list{3};
    if nargin < 4
        switch method
            case 'ratioStd'
                ratioTarget = 1/4;
            case 'localStd'
                increaseEps = 2;
        end
    else
        switch method
            case 'ratioStd'
                ratioTarget = str2double(arg_list{4});
            case 'localStd'
                increaseEps = str2double(arg_list{4});
        end
    end
end
switch method
    case 'ratioStd'
        fprintf(2,'Tonemapping: using method "%s" and ratioTarget = %f.\n',...
            method,ratioTarget)
    case 'localStd'
        fprintf(2,'Tonemapping: using method "%s" and increaseEps = %f.\n',...
            method,increaseEps)
end

%%% Add vpp mex files to path
addpath([fileparts(mfilename('fullpath')) '/../vpp-mex']);
addpath(fileparts(mfilename('fullpath')));

%%% Load required Octave packages
pkg load statistics
pkg load image

%%% Input
inputHandle = vpp_init_input(inputName);
if ~inputHandle, error('cannot load input %s',inputName); end

%%% Define size of output
H = inputHandle(1); % height of input
W = inputHandle(2); % width of input
D = 3;              % depth: output has color

%%% Output
outputHandle = vpp_init_output(outputName, [H,W,D]);
if ~outputHandle, error('cannot load output %s',outputName); end

%%% Get colormap (create variable "cmap")
warning ('off', 'Octave:data-file-in-path');
load([fileparts(mfilename('fullpath')) '/colormap.mat']);
if ~exist('cmap','var'), error('Couldn''t load the colormap.'); end

%%% MGF parameters
epsSqrt = 20;   % initialize epsilon. Its value is set automatically afterwards
radius  = 1;    % small images, need to avoid contrast halo
nIt     = 1;    % for speed
nScales = 0;    % filter all scales
gamma   = 1;    % evolution of gamma over the scales: eps <-- eps / scale^gamma
                %   - with gamma>0, reduce eps for low scales
                %   - with gamma=0, keep eps constant

%%% Get first frame; then loop until there is no frame available
u = double(vpp_read_frame(inputHandle));
n = 1;

%%% Initialiations
switch method
    case 'ratioStd'
        ratioThreshold = .05;
        step = NaN;

    case 'localStd'
        mask = u~=0;
        w = 4;
        k = ones(w*2+1) / (w*2+1)^2;
        localVar = imfilter(u.^2,k,'symmetric') - imfilter(u,k,'symmetric').^2;
        localStd = sqrt(max(0, localVar(mask) ));
        epsSqrt  = quantile(localStd,.75)*increaseEps;
end

while ~isscalar(u)

    %%% compute only statistics on valid pixels
    %%% u=0 happens due to the stabilization
    mask = u~=0;

    %%% Place input in a range where squared values does not explode
    %%% MGF does not need the input to be in a specific range, yet large values
    %%% cause numerical errors in single precision.
    %%% Use u = (u-mean(u(:)))/std(u(:)) if needed, or use double.
    m = mean(u(mask));
    u = u - m;

    %%% Replace stab borders by 0 to reduce creation of halos (u centered in 0)
    u(~mask) = 0;

    %%% Compute parameter epsilon, and base + detail decomposition
    switch method
        case 'ratioStd'
            %%% Estimation of epsilon: initializations
            uStd  = 1.4826*mad(u(mask),1); % = std(u(:));
            if uStd==0, uStd = eps; end % avoid ratio = NaN
            ratio = inf;
            nIter = 1;

            %%% Compute parameter epsilon, and base + detail decomposition
            while abs(ratio - ratioTarget) > ratioThreshold

                %%% Compute decomposition
                b = MGFoctave(u, u, epsSqrt^2, radius, nIt, nScales, gamma);
                d = u - b;

                %%% Compute ratio std(detail)/std(image)
                dStd  = 1.4826*mad(d(mask),1); % = std(d(:));
                ratio = dStd / uStd;

                % fprintf(2,...
                %     'Img #%04d\tIter #%d:\tEpsilon = %f\tRatio = %f\tStep = %f\n',...
                %     n, nIter, epsSqrt, ratio, step);

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
                if epsSqrt > 150
                    epsSqrt = 150;
                    if epsSqrtPrev == 150, break; end
                end

                %%% If estimation does not converge
                %%% This should not happen because of the condition above
                if nIter > 15, break; end
            end

        case 'localStd'
            %%% Local std (for valid pixels only)
            nIter = 2;
            localVar = imfilter(u.^2,k,'symmetric') - imfilter(u,k,'symmetric').^2;
            localStd = sqrt(max(0, localVar(mask) ));
            epsSqrt  = quantile(localStd,.75)*increaseEps;
            b = MGFoctave(u, u, epsSqrt^2, radius, nIt, nScales, gamma);
            d = u - b;

        otherwise
            error('Incorrect method.');
    end

    %%% Compress the base layer
    bMed = median(b(mask));
    bMad = median(abs(u(mask) - bMed)); % same as mad(b(mask),1);
    b    = max(-.5,min(+.5, (b-bMed) ./ (36 * 1.4826*bMad) )) / 2;

    %%% Compress the detail layer
    dMad = mad(d(mask),1);
    d    = d ./ (9 * 1.4826*dMad);

    %%% Gray tone-mapped image
    %%% We prefer clipping in black than in white, hence the .4 instead of .5
    v = .4 + b + d;

    %%% Count black and white clipping
    vClipBlack = sum( v(:) < 0 ) / numel(v);
    vClipWhite = sum( v(:) > 1 ) / numel(v);
    fprintf(2,...
        ['Img #%04d\tEps = %.1f (%d iter)' ...
         '\tClipped %.2f%% in black, %.2f%% in white.\n'],...
        n, epsSqrt, nIter-1, vClipBlack*100, vClipWhite*100);

    %%% Clip
    v = max(0,min(1, v ));
    v(~mask) = 0; % stabilization borders

    %%% Compute color coefficients using the (normalized+clipped) input image
    u = max(0,min(1, .5 + (u-bMed)./(9*1.4826*bMad) ));
    u = ind2rgb(uint16( u*(2^16-1) ), cmap);
    c = u ./ (sum(u,3)/3);

    %%% Apply color coefficient to vGray to get final colored tonemapped image
    v = 255 * v .* c;

    %%% Write frame
    if ~vpp_write_frame(outputHandle, v), break; end

    %%% Read next frame
    u = double(vpp_read_frame(inputHandle));
    n = n+1;
end

vpp_close(inputHandle);
vpp_close(outputHandle);
