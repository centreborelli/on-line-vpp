function R = MGFoctave( I, G, epsilon, radius, nIt, nScales, gamma )
% MGFoctave apply the multi-scale guided filter.
% R = MGFoctave( I, G, epsilon, radius, nIt, nScales, gamma )
%
% Charles Hessel -- CMLA, ENS-Paris-Saclay, 2018

if nScales==0
    nScales = max(1, floor(log2(min(size(I,1),size(I,2)))) );
elseif nScales < 0
    nScales = max(1, floor(log2(min(size(I,1),size(I,2)))) -nScales);
end
if nargin == 6
    gamma = 1;
end

if ( size(I,3) == 3 )
    error('The input image I must have only one channel.');
end

% switch between color and gray cases
color = ( size(G,3) == 3 );

% Use the filter defined in pyramid_filter (T. Mertens functions)
filter = pyramid_filter;

% Compute gaussian and laplacian pyramids of the input and guide images.
[~,I_lpyr] = glpyr(I,nScales);
G_gpyr = gpyr(G,nScales);

% Initialise the output R with the residual of the laplacian pyramid
R = I_lpyr{nScales};

% Filter the input image scale by scale starting from the coarsest one
% fprintf('Scale %02d: ',nScales);
if color
    R = igfcV2( R, G_gpyr{nScales}, epsilon/nScales^gamma, radius, nIt );
else
    R = igfV2(  R, G_gpyr{nScales}, epsilon/nScales^gamma, radius, nIt );
end

% Continue to recursively upsample, add coefficients and filter the image.
for scale = nScales-1:-1:1
    
    % Verify if we need to remove pixels after upsample
    odd = 2*size(R) - size(I_lpyr{scale});
    
    % Upsample previously filtered image, add Laplacian coefficients
    R = fastupsample(R,odd,filter) + I_lpyr{scale};
    
    % Filter
    % fprintf('Scale %02d: ',scale);
    if color
        R = igfcV2( R, G_gpyr{scale}, epsilon/scale^gamma, radius, nIt );
    else
        R = igfV2(  R, G_gpyr{scale}, epsilon/scale^gamma, radius, nIt );
    end
end

% This is a 1-dimensional 5-tap low pass filter. It is used as a 2D separable low
% pass filter for constructing Gaussian and Laplacian pyramids.
%
% tom.mertens@gmail.com, August 2007

function f = pyramid_filter
f = [.0625, .25, .375, .25, .0625];

% Simultaneous construction of Gaussian and Laplacian pyramid.
% Save computations when both are needed.
% Use gpyr to compute the gaussian pyramid only.
%
% This is base on a function by Tom Mertens.
%
% Charles Hessel, 2018

function [gpyr,lpyr] = glpyr(I,nlev)

% recursively build pyramid
gpyr = cell(nlev,1);
lpyr = cell(nlev,1);
filter = pyramid_filter;
J = I;
for l = 1:nlev - 1
    % apply low pass filter, and downsample
    I = fastdownsample(J,filter);
    gpyr{l} = J;
    % for each dimension, check if the upsampled version has to be odd
    odd = 2*size(I) - size(J);
    % in each level, store diff between image and upsampled low pass version
    lpyr{l} = J - fastupsample(I,odd,filter);
    J = I; % continue with low pass image
end
gpyr{nlev} = J; % the coarest level contains the residual low pass image
lpyr{nlev} = J;

function gpyr = gpyr(I,nlev)

% recursively build pyramid
gpyr = cell(nlev,1);
filter = pyramid_filter;
J = I;
for l = 1:nlev - 1
    % apply low pass filter, and downsample
    I = fastdownsample(J,filter);
    gpyr{l} = J;
    J = I; % continue with low pass image
end
gpyr{nlev} = J; % the coarest level contains the residual low pass image

% Upsampling procedure.
%
% Argments:
%   'I': greyscale image
%   'odd': 2-vector of binary values, indicates whether the upsampled image
%   should have odd size for the respective dimensions
%   'filter': upsampling filter
%
% If image width W is odd, then the resulting image will have width (W-1)/2+1,
% Same for height.
%
% tom.mertens@gmail.com, August 2007
% Charles Hessel, 2018

function R = fastupsample(I,odd,filter)

% increase resolution
I = padarray(I,[1 1 0],'replicate'); % pad the image with a 1-pixel border
r = 2*size(I,1);
c = 2*size(I,2);
k = size(I,3);

% horizontal, then vertical
Rt = zeros(r/2,c,k); Rt(:, 1:2:c, :) = 4*I; Rt = imfilter(Rt,filter);
R  = zeros(r  ,c,k);  R(1:2:r, :, :) = Rt;   R = imfilter(R,filter');

% remove the border
R = R(3:r - 2 - odd(1), 3:c - 2 - odd(2), :);

function R = fastdownsample(I, filter)

border_mode = 'symmetric';

% low pass, convolve with separable filter
R = imfilter(I,filter ,border_mode); %horizontal

% decimate in already filtered direction
R = R(:,1:2:end,:);

% low pass, convolve with separable filter
R = imfilter(R,filter',border_mode); %vertical

% decimate
R = R(1:2:end,:,:);

function R = igfcV2( I, G, epsilon, radius, T )
% igfcV2 is the iterated guided filter V2 with color guide
% (igfcV1 with unmodified guide).
%
% R = igfcV2( I, G, epsilon, radius, T );
% where T is the number of iterations

% fprintf('igfcV2: ');
% m = textprogressbar( 0 ,0, 50 );
R = I; % initialisation

for t=1:T
    % m = textprogressbar( t/T ,m, 50 );
    R = gfc( R, epsilon, radius, G );
end

function R = igfV2( I, G, epsilon, radius, T )

% igfV2 is the iterated guided filter V2 (igfV1 with unmodified guide).
% see thesis document. 
% R = igfV2( I, G, epsilon, radius, T );
% where T is the number of iterations
%
% \begin{equation}
%   R_{\text{IGF/v2}}^t = \bar{a}^t(p) I(p) + \bar{b}^t(p)
% \end{equation}
% where
% \begin{equation}
%   \bar{a}^t(p) = \mean \left\{
%     \frac{ \cov\{ I, R^{t-1} \} }{ \var\{ I \} + \epsilon } \right\}
% \end{equation}
% and
% \begin{equation}
%   \bar{b}^t = \mean \left\{ \mean\{ R^{t-1} \} - a^t(p) I(p)  \right\} ~,
% \end{equation}
% with $R_{\text{IGF/v2}}^0 = I$.

% fprintf('igfV2: ');
% m = textprogressbar( 0 ,0, 50 );
R = I; % initialisation

for t=1:T
    % m = textprogressbar( t/T ,m, 50 );
    R = gf( R, epsilon, radius, G );
end

function R = gf( I, epsilon, radius, G )
% gf is the guided filter.
% R = gf( I, epsilon, radius );
% OR
% R = gf( I, epsilon, radius, G ); for G \neq I.

%%% filer's shape for computation of the mean, var, covar.
H = ones(2*radius+1) ./ (2*radius+1)^2; 

switch nargin
    
    case 3 % G = I
        I_mean = imfilter(I   , H, 'symmetric');
        I_var  = imfilter(I.^2, H, 'symmetric') ...
                 - I_mean.^2;
        a      = I_var ./ ( I_var + epsilon );
        b      = ( 1 - a ) .* I_mean ;
        a_mean = imfilter( a , H, 'symmetric');
        b_mean = imfilter( b , H, 'symmetric');     
        R      = a_mean .* I + b_mean;

    case 4 % G \neq I
        I_mean = imfilter(I   , H, 'symmetric');
        G_mean = imfilter(G   , H, 'symmetric');
        G_var  = imfilter(G.^2, H, 'symmetric') - G_mean.^2;
        IG_cov = imfilter(I.*G, H, 'symmetric') - I_mean .* G_mean;
        a      = IG_cov ./ ( G_var + epsilon );
        b      = I_mean - a .* G_mean;
        a_mean = imfilter( a , H, 'symmetric');
        b_mean = imfilter( b , H, 'symmetric');     
        R      = a_mean .* G + b_mean;
                
    otherwise
        error('Incorrect number of arguments.')
end

function R = gfc( I, epsilon, radius, G )
% gfc is the guided filter with a color guide.
% R = gfc( I, epsilon, radius, G )

%%% filer's shape for computation of the mean, var, covar.
H = ones(2*radius+1) ./ (2*radius+1)^2; 

% mean vector and covariance matrix of color guide
Gm    = imfilter(G, H, 'symmetric');
Gv_rr = imfilter(G(:,:,1).^2, H, 'symmetric') - Gm(:,:,1).^2 + epsilon;
Gv_gg = imfilter(G(:,:,2).^2, H, 'symmetric') - Gm(:,:,2).^2 + epsilon;
Gv_bb = imfilter(G(:,:,3).^2, H, 'symmetric') - Gm(:,:,3).^2 + epsilon;
Gv_rg = imfilter(G(:,:,1).*G(:,:,2), H, 'symmetric') - Gm(:,:,1).*Gm(:,:,2);
Gv_rb = imfilter(G(:,:,1).*G(:,:,3), H, 'symmetric') - Gm(:,:,1).*Gm(:,:,3);
Gv_gb = imfilter(G(:,:,2).*G(:,:,3), H, 'symmetric') - Gm(:,:,2).*Gm(:,:,3);
detS  = +Gv_rr.*( Gv_gg.*Gv_bb - Gv_gb.*Gv_gb ) ...
        -Gv_rg.*( Gv_rg.*Gv_bb - Gv_rb.*Gv_gb ) ...
        +Gv_rb.*( Gv_rg.*Gv_gb - Gv_rb.*Gv_gg ) ;

% mean of gray input and vector of covariance between input and guide
Im    = imfilter(I, H, 'symmetric');
IGc_r = imfilter(I.*G(:,:,1), H, 'symmetric') - Im.*Gm(:,:,1);
IGc_g = imfilter(I.*G(:,:,2), H, 'symmetric') - Im.*Gm(:,:,2);
IGc_b = imfilter(I.*G(:,:,3), H, 'symmetric') - Im.*Gm(:,:,3);

% solve the system and get vector of linear coefficients "a"
detR  = +IGc_r.*( Gv_gg.*Gv_bb - Gv_gb.*Gv_gb ) ...
        -Gv_rg.*( IGc_g.*Gv_bb - IGc_b.*Gv_gb ) ...
        +Gv_rb.*( IGc_g.*Gv_gb - IGc_b.*Gv_gg ) ;
detG  = +Gv_rr.*( IGc_g.*Gv_bb - IGc_b.*Gv_gb ) ...
        -IGc_r.*( Gv_rg.*Gv_bb - Gv_rb.*Gv_gb ) ...
        +Gv_rb.*( Gv_rg.*IGc_b - Gv_rb.*IGc_g ) ;
detB  = +Gv_rr.*( Gv_gg.*IGc_b - Gv_gb.*IGc_g ) ...
        -Gv_rg.*( Gv_rg.*IGc_b - Gv_rb.*IGc_g ) ...
        +IGc_r.*( Gv_rg.*Gv_gb - Gv_rb.*Gv_gg ) ;

a     = cat(3, detR./detS, detG./detS, detB./detS);

% average linear coefficients of overlapping windows
am    = imfilter(a, H, 'symmetric');
bm    = imfilter(Im - a(:,:,1).*Gm(:,:,1) ...
                    - a(:,:,2).*Gm(:,:,2) ...
                    - a(:,:,3).*Gm(:,:,3), H, 'symmetric');

% compute output gray image
R     =   am(:,:,1).*G(:,:,1)...
        + am(:,:,2).*G(:,:,2)...
        + am(:,:,3).*G(:,:,3)...
        + bm;

% function msg_length = textprogressbar( x ,prev_msg_length, barlength )
% 
% ers         = repmat('\b',1,prev_msg_length);
% n           = round(x*barlength);
% bar         = sprintf('[%s%s]\n',repmat('#',1,n),...
%                                  repmat(' ',1,barlength-n));
% msg_length  = length(bar);
% fprintf([ers bar]);
