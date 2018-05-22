function ops0 = init_ops(ops0)
%%
ops0.toolbox_path = 'C:\home\matlab_svn\Suite2P';
if exist(ops0.toolbox_path, 'dir')
    addpath(genpath(ops0.toolbox_path)) % add local path to the toolbox
else
    error('toolbox_path does not exist, please change toolbox_path');
end

% mex -largeArrayDims SpikeDetection/deconvL0.c (or .cpp) % MAKE SURE YOU COMPILE THIS FIRST FOR DECONVOLUTION

ops0.useGPU                 = 0; % if you can use an Nvidia GPU in matlab this accelerates registration approx 3 times. You only need the Nvidia drivers installed (not CUDA).

% root paths for files and temporary storage (ideally an SSD drive. my SSD is C:/)
ops0.RootStorage            =getOr(ops0, {'RootStorage'}, 'F:\home\ImagingData\DualLickMice\'); % Suite2P assumes a folder structure, check out README file
% ops0.RootStorage            = 'G:\Kosuke\Data\Fh\DualLickMice\'; % Suite2P assumes a folder structure, check out README file
ops0.temp_tiff              = 'G:\Kosuke\Data\temp.tif'; % copy each remote tiff locally first, into this file
ops0.RegFileRoot            = 'G:\Kosuke\Data';  % location for binary file, better to be SSD drive.
if ~exist(ops0.RegFileRoot,'dir')
    try mkdir(ops0.RegFileRoot)
    catch 
        ops0.RegFileRoot            = 'C:\tmp\';  % 
        if ~exist(ops0.RegFileRoot,'dir')
            mkdir(ops0.RegFileRoot);
            ops0.temp_tiff              = fullfile(ops0.RegFileRoot,'temp.tif'); %'G:\Kosuke\Data\temp.tif';
        end
    end
end
ops0.DeleteBin              = 1; % set to 1 for batch processing on a limited hard drive
ops0.ResultsSavePath        = ops0.RootStorage ; % a folder structure is created inside
ops0.RegFileTiffLocation    = ops0.RootStorage ; %'D:/DATA/'; % leave empty to NOT save registered tiffs (slow)

% registration options
% ops0.doRegistration         = 1; % skip (0) if data is already registered
ops0.doRegistration         = getOr(ops0, {'doRegistration'}, 1);
ops0.showTargetRegistration = 1; % shows the image targets for all planes to be registered
ops0.PhaseCorrelation       = 1; % set to 0 for non-whitened cross-correlation
ops0.SubPixel               = Inf; % 2 is alignment by 0.5 pixel, Inf is the exact number from phase correlation
ops0.NimgFirstRegistration  = 500; % number of images to include in the first registration pass 
ops0.nimgbegend             = 250; % frames to average at beginning and end of blocks
ops0.MaxMovementPixel       = 100; % If the estimated movement shift is more than this value, use interpolation. 
ops0.PhaseCorrBlurSTD       = 1.5; % To stabilize the selection of movement shift, filter the phase correlation with Gaussian filter.

% cell detection options
% ops0.clustModel             = 'neuropil'; % standard or neuropil, or fast_spectral (by KH)
ops0.clustModel             = 'HDBCellScan'; % The fastest hierachical density based cell scan approach (using python internally), by KH 20171102.
ops0.neuropilSub            = 'surround'; % none, surround or model
ops0.ShowCellMap            = 1; % during optimization, show a figure of the clusters
ops0.Nk0                    = 100; % how many clusters to start with
ops0.Nk                     = 26;  % how many clusters to end with (before anatomical segmentation)
ops0.sig                    = 0.5;  % spatial smoothing length in pixels; encourages localized clusters

% -- how many (binned) timepoints to do the SVD based on. When the total frame is 30,000, and NavgFramesSVD = 5000,
% it will calculate the average of each (30000/5000 =) 6 frames and use for SVD calculation.
ops0.NavgFramesSVD          = 3000; % 
ops0.nSVDMax                = 1000;  % max number of SVD components for cell clustering
ops0.getSVDcomps                    = getOr(ops0, {'getSVDcomps'}, 0);   % whether to use old SVD comp
ops0.nSVD                           = getOr(ops0, {'nSVD'}, 1000);   % how many SVD components to save to disk
ops0.nSVDforROI             = 100; % how many SVD components for cell clustering


% ops0.nFramesAvgForSVD       = 6;   % nFramesAvgForSVD = TotalFrames/NavgFramesSVD.
% if nFramesAvgForSVD is defined, it overwrites NavgFramesSVD and nSVDforROI
% NavgFramesSVD= TotalFrames/nFramesAvgForSVD,
% nSVDforROI   = min(1000,NavgFramesSVD);

% red channel options
% redratio = red pixels inside / red pixels outside
% redcell = redratio > mean(redratio) + redthres*std(redratio)
% notred = redratio < mean(redratio) + redmax*std(redratio)
ops0.redthres               = 1.5; % the higher the thres the less red cells
ops0.redmax                 = 1; % the higher the max the more NON-red cells

% spike deconvolution options
ops0.imageRate              = 30;   % imaging rate (cumulative over planes!). Approximate, for initialization of deconvolution kernel.
ops0.sensorTau              = 1; % decay half-life (or timescale). Approximate, for initialization of deconvolution kernel.
ops0.maxNeurop              = Inf; % for the neuropil contamination to be less than this (sometimes good, i.e. for interneurons)
ops0.recomputeKernel        = 1; % whether to re-estimate kernel during optimization (default kernel is "reasonable", if you give good timescales)
ops0.sameKernel             = 1; % whether the same kernel should be estimated for all neurons (robust, only set to 0 if SNR is high and recordings are long)


%% 
% ops0.nimgbegend                     = getOr(ops0, {'nimgbegend'}, 0);
ops0.splitROIs                      = getOr(ops0, {'splitROIs'}, 1);
ops0.LoadRegMean                    = getOr(ops0, {'LoadRegMean'}, 0);
ops0.NiterPrealign                  = getOr(ops0, {'NiterPrealign'}, 10);
ops0.registrationUpsample           = getOr(ops0, {'registrationUpsample'}, 1);  % upsampling factor during registration, 1 for no upsampling is much faster, 2 may give better subpixel accuracy
ops0.niterclustering                = getOr(ops0, {'niterclustering'}, 50);   % how many iterations of clustering
ops0.getROIs                        = getOr(ops0, {'getROIs'}, 1);   % whether to calculate ROI and save SVD components

% ops0.diameter                       = clustrules.diameter;

% ops0.mimg = mean(single(squeeze(IMG(yFOVs(:,jj),xFOVs(:,jj),ops.planesToProcess(ii),:))),3);
% use ops.dsprealign to align other images
% ops0.RefImg=ops0.mimg;
% ops0.IsRefChannel= 1;

% may not be the best solution, but for back-compativility.
ops0.PlaneID = getOr(ops0,'PlaneID',1);
ops0.ViewID  = getOr(ops0,'ViewID',1);
ops0.ChannelID  = getOr(ops0,'ChannelID',1);

% ops0.Ly   =  getOr(ops0,'Ly',512);
% ops0.Lx   =  getOr(ops0,'Lx',512);
