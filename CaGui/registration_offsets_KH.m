function [dv, corr, regdata] = registration_offsets_KH(data, ops, removeMean,varargin)
% updated 20171006: to use RefImage, use 4th argument. 
%% Parameters
[ly, lx, nFrames] = size(data);
if nargin>=4
    refImg = varargin{1};
else
    refImg = ops.mimg;
end

subpixel = getOr(ops, {'subPixel' 'SubPixel'}, 1); % subpixel factor
usFac = getOr(ops, 'registrationUpsample', 1); % factor to upsample
phaseCorrelation = getOr(ops, {'phaseCorrelation' 'PhaseCorrelation'}, false);
useGPU = getOr(ops, 'useGPU', false);

maskSlope   = getOr(ops, 'maskSlope', 2); % slope on taper mask preapplied to image. was 2, then 1.2
% SD pixels of gaussian smoothing applied to correlation map (MOM likes .6)
defaultSig = 1.15;
smoothSigma   = getOr(ops, 'sig', defaultSig); % default value is 1.15;
if smoothSigma ~= defaultSig, fprintf('SmoothSigma=%f\n',smoothSigma); end; % to notify user that this program is using non-default value
smoothSigma = smoothSigma/sqrt(usFac);

if nargout > 2 % translation required
    translate = true;
    fy = ifftshift((-fix(ly/2):ceil(ly/2) - 1)/ly)';% freq along first dimension
    fx = ifftshift((-fix(lx/2):ceil(lx/2) - 1)/lx); % freq along second dimension
else
    translate = false;
end
%% Prepare common arrays
lyus = usFac*ly;
lxus = usFac*lx;
% Taper mask
[ys, xs] = ndgrid(1:ly, 1:lx);
ys = abs(ys - mean(ys(:)));
xs = abs(xs - mean(xs(:)));
mY      = max(ys(:)) - 4;
mX      = max(xs(:)) - 4;
maskMul = single(1./(1 + exp((ys - mY)/maskSlope)) ./(1 + exp((xs - mX)/maskSlope)));
maskOffset = mean(refImg(:))*(1 - maskMul);

% fixed gaussian probability mask to limit movement correction: 
% introduced by KH 20190617
rs2 = (ys-(1+ly)/2).^2 + (xs-(1+lx)/2).^2;
gmaskVarDefault = inf;
gmaskVar   = getOr(ops, 'gmaskVar', gmaskVarDefault); % reduce the probability of picking up the movements correction outside the 30 pixels (~30um for Olympus x25 lens)
if gmaskVar ~=gmaskVarDefault,    fprintf('gmaskVar=%3.1f\n',gmaskVar); end
gmask_corr = single(exp(-rs2/gmaskVar)); % this mask is shifted so that it can be applied to the ffted correlation map. 

% Array indices for centre of mass clip window
[yClipRef, xClipRef] = ndgrid(-2:2, -2:2);
xClipRef = xClipRef(:);
yClipRef = yClipRef(:);
nClipPixels = numel(xClipRef);
% Array indices for embedding fourier components in a larger array
yEmbedRef = [1:fix((ly + 1)/2) (lyus - fix(ly/2) + 1):lyus];
xEmbedRef = [1:fix((lx + 1)/2) (lxus - fix(lx/2) + 1):lxus];
% Array indices for correlation clip window. Assumes at jitter +/-lCorr
lCorr = 50;
xCorrRef = [(usFac*lx - lCorr + 1):usFac*lx 1:(lCorr + 1)];
yCorrRef = [(usFac*ly - lCorr + 1):usFac*ly 1:(lCorr + 1)];
% Smoothing filter in frequency domain
hgx = exp(-(((0:lx-1) - fix(lx/2))/smoothSigma).^2);
hgy = exp(-(((0:ly-1) - fix(ly/2))/smoothSigma).^2);
hg = hgy'*hgx;
fhg = real(fftn(ifftshift(single(hg/sum(hg(:))))));
% Prepare data arrays
cfRefImg = conj(fftn(refImg));
eps0 = single(1e-20);
if phaseCorrelation
    cfRefImg = cfRefImg./(eps0 + abs(cfRefImg)).*fhg;
end

if useGPU
    batchSize = getBatchSize(lyus*lxus);
    maskMul = gpuArray(maskMul);
    maskOffset = gpuArray(maskOffset);
    cfRefImg = gpuArray(cfRefImg);
    eps0 = gpuArray(eps0);
    corrUps = zeros(lyus, lxus, batchSize, 'single', 'gpuArray');
    if nargout > 2
        fx = gpuArray(fx);
        fy = gpuArray(fy);
    end
else
    batchSize = 100;
    corrUps = zeros(lyus, lxus, batchSize, 'single');
end
%% Work through data in batches
dv = zeros(nFrames, 2);
corr = zeros(nFrames, 1);
if translate
    regdata = zeros(ly, lx, nFrames, 'single');
end
nBatches = ceil(nFrames/batchSize);
for bi = 1:nBatches
    fi = (bi - 1)*batchSize + 1:min(bi*batchSize, nFrames);
    if bi == nBatches
        % the last batch will usually have less frames
        corrUps = corrUps(:,:,1:numel(fi));
    end
    if useGPU
        batchData = gpuArray(single(data(:,:,fi)));
    else
        batchData = single(data(:,:,fi));
    end
    corrMap = fft2(bsxfun(@plus, maskOffset, bsxfun(@times, maskMul, batchData)));
    if phaseCorrelation
        corrMap = bsxfun(@times, corrMap./(eps0 + abs(corrMap)), cfRefImg);
    else
        corrMap = bsxfun(@times, corrMap, cfRefImg);
    end
    % embed in a larger array and compute 2D inverse fft to get correlation map
    corrUps(yEmbedRef,xEmbedRef,:) = corrMap;
    corrUps = real(ifft2(corrUps));
    % to reduce the probability of moving too far (added by KH 20190617)
    corrUps=bsxfun(@times,corrUps,gmask_corr);
    corrClip = corrUps(yCorrRef,xCorrRef,:);
    
    % added by Marius 20.07.2016, smooth the correlation maps
    %     corrClipSmooth = my_conv2(corrClip, 1, [1 2]);
    
    % added by Kosuke 23.11.2016 to make it option.
    DefaultPhaseCorrBlurSTD=1.5;
    ConvSTD = getOr(ops, {'PhaseCorrBlurSTD'}, DefaultPhaseCorrBlurSTD);
    
    if ConvSTD~=DefaultPhaseCorrBlurSTD, fprintf('%s=%d','PhaseCorrBlurSTD',ConvSTD); end
    corrClipSmooth = my_conv2(corrClip, ConvSTD, [1 2]);
%     corrClipSmooth = bsxfun(@times,corrClipSmooth,gmask_corr);
  
    % find peak
    [dmax, iy] = max(corrClipSmooth, [], 1);
    iy = gather_try(iy);
    dmax = gather_try(dmax);
    [dmax, ix] = max(dmax, [], 2);
    iy = reshape(...
        iy(sub2ind([size(iy,2) size(iy,3)], ix(:), (1:size(iy,3))')),...
        1, 1, []);
    if subpixel > 1
        iy = min(max(iy, 3), 2*lCorr - 1);
        ix = min(max(ix, 3), 2*lCorr - 1);
        clipX = bsxfun(@plus, xClipRef', ix);
        clipY = bsxfun(@plus, yClipRef, iy);
        clipF = reshape(repmat(1:size(clipX, 3), nClipPixels, 1), [], 1);
        cczoom = reshape(...
            gather_try(corrClip(sub2ind(size(corrClip), clipY(:), clipX(:), clipF))),...
            nClipPixels, 1, []);
        bcorr = sum(cczoom, 1);
        cczoom = bsxfun(@rdivide, cczoom, bcorr);
        ix = ix + sum(bsxfun(@times, xClipRef, cczoom), 1);
        iy = iy + sum(bsxfun(@times, yClipRef, cczoom), 1);
    else
        bcorr = dmax;
    end
    ix = (ix - lCorr - 1)/usFac;
    iy = (iy - lCorr - 1)/usFac;
    if isfinite(subpixel)
        ix = round(subpixel*ix)./subpixel;
        iy = round(subpixel*iy)./subpixel;
    end
    if translate % do translation using registration offsets in fourier domain
        phaseShift = bsxfun(@times,...
            exp(1j*2*pi*bsxfun(@times, fy, iy)),... y rotation
            exp(1j*2*pi*bsxfun(@times, fx, ix))); % x rotation
        res = real(ifft2(fft2(batchData).*phaseShift));
        regdata(:,:,fi) = gather_try(res);
    end
    dv(fi,:) = [iy(:) ix(:)];
    corr(fi) = squeeze(bcorr);
end
%% Post-processing
if nargin > 2 && removeMean
    dv = bsxfun(@minus, dv, mean(dv,1));
end
