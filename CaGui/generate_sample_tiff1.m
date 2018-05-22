function F = generate_sample_tiff1()
% generate 2p-sample data
M = 512;
N = 512;
Fs = 33;
Tmax = 7;
T = linspace(0,Tmax,Tmax*Fs);
Tlen = length(T);

x = [1:M];
y = [1:N];
[xgrid,ygrid]=meshgrid(x,y);
CaTau = 0.6*Fs;

%% cell data
ii=1;
n(ii).x = round(M/4); % x position of a cell
n(ii).y = round(N/4); % y position of a cell
n(ii).r = 12; % radius 
n(ii).mask = exp(-((xgrid-n(ii).x).^2+(ygrid-n(ii).y).^2)/(2*n(ii).r.^2));
n(ii).mask = n(ii).mask .* (n(ii).mask>0.5);

n(ii).spk = 30*sparse(rand(1,Tlen)<1/Fs);
n(ii).F  = 0;


ii=ii+1;
n(ii).x = round(M*0.75); % x position of a cell
n(ii).y = round(N*0.25); % y position of a cell
n(ii).r = 15; % radius 
n(ii).mask = exp(-((xgrid-n(ii).x).^2+(ygrid-n(ii).y).^2)/(2*n(ii).r.^2));
n(ii).mask = n(ii).mask .* (n(ii).mask>0.5);

n(ii).spk = 30*sparse(rand(1,Tlen)<2/Fs);
n(ii).F  = 0;



ii=ii+1;
n(ii).x = round(M*0.25); % x position of a cell
n(ii).y = round(N*0.75); % y position of a cell
n(ii).r = 18; % radius 
n(ii).mask = exp(-((xgrid-n(ii).x).^2+(ygrid-n(ii).y).^2)/(2*n(ii).r.^2));
n(ii).mask = n(ii).mask .* (n(ii).mask>0.5);

n(ii).spk = 30*sparse(rand(1,Tlen)<4/Fs);
n(ii).F  = 0;



ii=ii+1;
n(ii).x = round(M*0.75); % x position of a cell
n(ii).y = round(N*0.75); % y position of a cell
n(ii).r = 15; % radius 
n(ii).mask = exp(-((xgrid-n(ii).x).^2+(ygrid-n(ii).y).^2)/(2*n(ii).r.^2));
n(ii).mask = n(ii).mask .* (n(ii).mask>0.5);

n(ii).spk = 30*sparse(rand(1,Tlen)<0.5/Fs);
n(ii).F  = 0;


F = zeros(N,M,Tlen);
F = uint16(F);

BaseNoiseSigma = 3;
%% movement artifact
theta = 40/180*pi;
movement_vector = [cos(theta), sin(theta)];
movement_pulse = 10*sparse(rand(1,Tlen)<0.3/Fs);
movement_trace = 0;
random_movement = 2;
%%
myfigure('2P simulated data');clf;
imgh=image(F(:,:,1));
colormap gray;

ds = zeros(2,Tlen);
tmp_old = F(:,:,1);

for tt=1:Tlen
    tmp = uint16(zeros(N,M));
    for nn=1:length(n)
        if tt==1
            OldF = 0;
        else
            OldF = n(nn).F(tt-1);
        end
        n(nn).F(tt) = (1-1/CaTau)*OldF + n(nn).spk(tt); % CaTau * dF/dT = -F + spk
        movement_trace = (1-1/(0.3*Fs))*movement_trace + movement_pulse(tt);
        tmp = tmp+uint16(n(nn).mask .* n(nn).F(tt) .* rand(N,M)*2);
%         F(:,:,tt) =  F(:,:,tt) + uint16();
    end
    
    tmp= tmp+ uint16(rand(N,M)*BaseNoiseSigma + BaseNoiseSigma*10*sprand(N,M,0.01));
    
    dx = round(movement_trace*movement_vector(1)+random_movement*randn(1));
    dy = round(movement_trace*movement_vector(2)+random_movement*randn(1));
    F(:,:,tt) = circshift(tmp,[dy,dx]);
    ds(:,tt)=[dx,dy]';
%     tmp_old = tmp;
    if rem(tt,10)==0
        set(imgh,'CData',F(:,:,tt)); 
        title(sprintf('t=%3.1f',T(tt)));
        drawnow;
    end
end

%% 
% F is 16bit integer image stacks.
% FilePath = 'G:\Kosuke\Data\F\SimulatedData\20161206\1';
FilePath = 'D:\home\ImagingData\DualLickMice\SimulatedData\20161206\1';
mkdir(FilePath);
if isempty(mfilename)
    TiffFileName = fullfile(FilePath,['generate_sample_tiff1.tif']);
else  
    TiffFileName = fullfile(FilePath,[mfilename,'.tif']);
end

[FilePath,FileName,Ext]=fileparts(TiffFileName);

fprintf('Writing Tiff data in %s...\n',TiffFileName)
bitspersamp = 16;
 TiffWriter(F,TiffFileName,bitspersamp);
 save(fullfile(FilePath,[FileName,'.mat']),'F','ds','n');
 
fprintf('Writing matlab data in %s...\n',fullfile(FilePath,[FileName,'.mat']));
 
 %% run  master_file_example_kh 
 
 %% SET ALL DEFAULT OPTIONS HERE
% check out the README file for detailed instructions (and extra options)
addpath('F:\home\ImagingData\DualLickMice\') % add the path to your make_db file
% addpath('D:\home\ImagingData\DualLickMice'); % For the case of using external HD. Change Drive letter as you need.

db = [];
% overwrite any of these default options in your make_db file for individual experiments
% make_db_example; % RUN YOUR OWN MAKE_DB SCRIPT TO RUN HERE
% make_db_B6J448;
% make_db_B6J449;
make_db_SimulatedData;

ops0.toolbox_path = 'C:\home\matlab_svn\Suite2P';
if exist(ops0.toolbox_path, 'dir')
	addpath(genpath(ops0.toolbox_path)) % add local path to the toolbox
else
	error('toolbox_path does not exist, please change toolbox_path');
end

% mex -largeArrayDims SpikeDetection/deconvL0.c (or .cpp) % MAKE SURE YOU COMPILE THIS FIRST FOR DECONVOLUTION

ops0.useGPU                 = 1; % if you can use an Nvidia GPU in matlab this accelerates registration approx 3 times. You only need the Nvidia drivers installed (not CUDA).

% PCEnv = 'mynotePC';
PCEnv = 'golgi';

switch PCEnv
    case 'golgi'
        % root paths for files and temporary storage (ideally an SSD drive. my SSD is C:/)
        ops0.RootStorage            = 'F:\home\ImagingData\DualLickMice\'; % Suite2P assumes a folder structure, check out README file
        ops0.temp_tiff              = 'G:\Kosuke\Data\temp.tif'; % copy each remote tiff locally first, into this file
        ops0.RegFileRoot            = 'G:\Kosuke\Data';  % location for binary file, better to be SSD drive.
        ops0.DeleteBin              = 1; % set to 1 for batch processing on a limited hard drive
        ops0.ResultsSavePath        = 'G:\Kosuke\Data\F'; % a folder structure is created inside
        ops0.RegFileTiffLocation    = 'G:\Kosuke\Data\F'; %'D:/DATA/'; % leave empty to NOT save registered tiffs (slow)
        
    case 'mynotePC'
        % root paths for files and temporary storage (ideally an SSD drive. my SSD is C:/)
        ops0.RootStorage            = 'D:\home\ImagingData\DualLickMice\'; % Suite2P assumes a folder structure, check out README file
        ops0.temp_tiff              = 'D:\home\ImagingData\temp.tif'; % copy each remote tiff locally first, into this file
        ops0.RegFileRoot            = 'D:\home\ImagingData\';  % location for binary file, better to be SSD drive.
        ops0.DeleteBin              = 1; % set to 1 for batch processing on a limited hard drive
        ops0.ResultsSavePath        = 'D:\home\ImagingData\'; % a folder structure is created inside
        ops0.RegFileTiffLocation    = 'D:\home\ImagingData\'; %'D:/DATA/'; % leave empty to NOT save registered tiffs (slow)
end
% registration options
ops0.doRegistration         = 1; % skip (0) if data is already registered
ops0.showTargetRegistration = 1; % shows the image targets for all planes to be registered
ops0.PhaseCorrelation       = 1; % set to 0 for non-whitened cross-correlation
ops0.SubPixel               = Inf; % 2 is alignment by 0.5 pixel, Inf is the exact number from phase correlation
ops0.NimgFirstRegistration  = 50; % number of images to include in the first registration pass 
ops0.nimgbegend             = 150; % frames to average at beginning and end of blocks
ops0.PhaseCorrBlurSTD       = 1.2; % To stabilize the selection of movement shift, filter the phase correlation with Gaussian filter.
ops0.MaxMovementPixel       = 80; % If the estimated movement shift is more than this value, use interpolation. 

% cell detection options
ops0.clustModel             = 'neuropil'; % standard or neuropil
ops0.neuropilSub            = 'model'; % none, surround or model
ops0.ShowCellMap            = 1; % during optimization, show a figure of the clusters
ops0.Nk0                    = 8; % how many clusters to start with
ops0.Nk                     = 4;  % how many clusters to end with (before anatomical segmentation)
ops0.sig                    = 0.5;  % spatial smoothing length in pixels; encourages localized clusters
ops0.nSVDforROI             = 80; % how many SVD components for cell clustering
ops0.NavgFramesSVD          = 160; % how many (binned) timepoints to do the SVD based on
clustrules.diameter         = 40; % expected diameter of cells (used for 0.25 * pi/4*diam^2 < npixels < 10*pi/4*diam^2)

% red channel options
% redratio = red pixels inside / red pixels outside
% redcell = redratio > mean(redratio) + redthres*std(redratio)
% notred = redratio < mean(redratio) + redmax*std(redratio)
ops0.redthres               = 1.5; % the higher the thres the less red cells
ops0.redmax                 = 1; % the higher the max the more NON-red cells

% spike deconvolution options
ops0.imageRate              = 33;   % imaging rate (cumulative over planes!). Approximate, for initialization of deconvolution kernel.
ops0.sensorTau              = 2; % decay half-life (or timescale). Approximate, for initialization of deconvolution kernel.
ops0.maxNeurop              = Inf; % for the neuropil contamination to be less than this (sometimes good, ii.e. for interneurons)
ops0.recomputeKernel        = 1; % whether to re-estimate kernel during optimization (default kernel is "reasonable", if you give good timescales)
ops0.sameKernel             = 1; % whether the same kernel should be estimated for all neurons (robust, only set to 0 if SNR is high and recordings are long)

db0 = db;
%%  this is the half of run_pipeline(db(iexp), ops0, clustrules);

ops0.nimgbegend                     = getOr(ops0, {'nimgbegend'}, 0);
ops0.splitROIs                      = getOr(ops0, {'splitROIs'}, 1);
ops0.LoadRegMean                    = getOr(ops0, {'LoadRegMean'}, 0);
ops0.NiterPrealign                  = getOr(ops0, {'NiterPrealign'}, 10);
ops0.registrationUpsample           = getOr(ops0, {'registrationUpsample'}, 1);  % upsampling factor during registration, 1 for no upsampling is much faster, 2 may give better subpixel accuracy
ops0.niterclustering                = getOr(ops0, {'niterclustering'}, 50);   % how many iterations of clustering
ops0.getROIs                        = getOr(ops0, {'getROIs'}, 1);   % whether to run the optimization
ops0.getSVDcomps                    = getOr(ops0, {'getSVDcomps'}, 0);   % whether to save SVD components to disk for later processing
ops0.nSVD                           = getOr(ops0, {'nSVD'}, 1000);   % how many SVD components to save to disk

ops0.diameter                       = clustrules.diameter;

clustrules = get_clustrules(clustrules);
 
ops = build_ops3(db(end), ops0);


clustModel     = getOr(ops, {'clustModel'}, 'standard');
neuropilSub    = getOr(ops, {'neuropilSub'}, 'surround');
splitBlocks    = getOr(ops, {'splitBlocks'}, 'none');
 
   ops1         = reg2P(ops);  % do registration
   
%% confirm how good the movement correction is. 
myfigure('Movement Correction Result')
ds_estim=   ops1{1}.DS;
clf;
subh(1)=subplot(3,1,1);
plot(ds(1,:),'k'); hold on;
plot(ds_estim(:,2),'r--');

subh(2)=subplot(3,1,2);
plot(ds(2,:),'k'); hold on;
plot(ds_estim(:,1),'r--');

subh(3)=subplot(3,1,3);
plot(ops1{1}.CorrFrame);
plot(sqrt(sum(ops1{1}.DS.^2,2)));
linkaxes(subh,'x')
%% first, obtain U which is U*V = F. 
%  U is eigenvalues which represents 
ii=1;
iplane  = ops.planesToProcess(ii);
ops     = ops1{ii};

ops.iplane  = iplane;

if getOr(ops, {'getSVDcomps'}, 0)
    ops    = get_svdcomps(ops);
end

if ops.getROIs || getOr(ops, {'writeSVDroi'}, 0)
    [ops, U, Sv]    = get_svdForROI(ops);
end

myfigure('Obtained EigenImage');clf;
NN=3;
for ii=1:NN^2
    mysubplot(NN,NN,ii); imagesc(U(:,:,ii));
    box off;
    axis off;
end 
%% 
 if ops.getROIs
        switch clustModel
            case 'standard'
                [ops, stat, res]  = fast_clustering_kh01(ops,U, Sv);
            case 'neuropil'                    
%                 [ops, stat, res]  = fast_clustering_with_neuropil(ops,U, Sv);
                  % better model of the neuropil
                  [ops, stat, res]  = fastClustNeuropilCoef(ops,U, Sv);
        end
                
        [stat2, res2] = apply_ROIrules(ops, stat, res, clustrules);
        
        switch neuropilSub
            case 'surround'
                get_signals_and_neuropil(ops, iplane);
            case 'none'
                get_signals(ops, iplane);
            case 'model'
                get_signals_NEUmodel(ops, iplane);
        end
    end
    

    if ops.DeleteBin
        fclose('all');
        delete(ops.RegFile);        % delete temporary bin file
    end


%%
LoadName = sprintf('%s/F_%s_%s_plane%d_Nk%d.mat',ops.ResultsSavePath,ops.mouse_name, ops.date, iplane, ops.Nk);
load(LoadName,'ops', 'res', 'stat', ...
        'stat0', 'res0', 'Fcell', 'FcellNeu', 'clustrules');

%%    

myfigure('Detected Cell');clf;
subplot(2,1,1);
imagesc(ops.mimg1);

subplot(2,1,2);
cellind = find([stat(:).igood]);
F = Fcell{1}(cellind,:)';
zF = zscore(F,0,1);
plot(zF,'LineWidth',2);

trueF = reshape([n(:).F],length(n(1).F),length(n));
z_trueF = zscore(trueF,0,1);
% for ii=cellind
hold on;
plot(z_trueF,'k--'); 
% end
[SaveDir,FileName,Ext]=fileparts(TiffFileName);

print(gcf,'-dpng',fullfile(SaveDir,[FileName,'.png']));
print(gcf,'-depsc',fullfile(SaveDir,[FileName,'.eps']));
