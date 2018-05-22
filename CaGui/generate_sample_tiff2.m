function F = generate_sample_tiff2()
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
% FilePath = 'D:\home\ImagingData\DualLickMice\SimulatedData\20161206\1';
FilePath = 'C:\home\ImagingExperiments\SimulatedData\SimulatedOne\20171020';
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
% addpath('F:\home\ImagingData\DualLickMice\S') % add the path to your make_db file
% addpath('D:\home\ImagingData\DualLickMice'); % For the case of using external HD. Change Drive letter as you need.
addpath('C:\home\ImagingExperiments\SimulatedData'); % for T470sKH 

db = [];
% overwrite any of these default options in your make_db file for individual experiments
% make_db_example; % RUN YOUR OWN MAKE_DB SCRIPT TO RUN HERE
% make_db_B6J448;
% make_db_B6J449;
[db,ops0]=make_db_SimulatedData;

ops0.toolbox_path = 'C:\home\matlab_svn\Suite2P';
if exist(ops0.toolbox_path, 'dir')
	addpath(genpath(ops0.toolbox_path)) % add local path to the toolbox
else
	error('toolbox_path does not exist, please change toolbox_path');
end

% mex -largeArrayDims SpikeDetection/deconvL0.c (or .cpp) % MAKE SURE YOU COMPILE THIS FIRST FOR DECONVOLUTION

ops0.useGPU                 = 0; % if you can use an Nvidia GPU in matlab this accelerates registration approx 3 times. You only need the Nvidia drivers installed (not CUDA).

PCEnv = 'mynotePC';
% PCEnv = 'golgi';

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
%         C:\home\ImagingExperiments\SimulatedData\SimulatedOne\20171020\1
        % root paths for files and temporary storage (ideally an SSD drive. my SSD is C:/)
        ops0.RootStorage            = 'C:\home\ImagingExperiments\SimulatedData\'; % Suite2P assumes a folder structure, check out README file
        ops0.temp_tiff              = 'C:\home\ImagingExperiments\temp.tif'; % copy each remote tiff locally first, into this file
        ops0.RegFileRoot            = 'C:\home\ImagingExperiments\';  % location for binary file, better to be SSD drive.
        ops0.DeleteBin              = 1; % set to 1 for batch processing on a limited hard drive
        ops0.ResultsSavePath        = 'C:\home\ImagingExperiments\SimulatedData\'; % a folder structure is created inside
        ops0.RegFileTiffLocation    = 'C:\home\ImagingExperiments\SimulatedData\'; %'D:/DATA/'; % leave empty to NOT save registered tiffs (slow)
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
ops0.writeSVDroi                    = getOr(ops0, {'writeSVDroi'}, 1);
ops0.diameter                       = clustrules.diameter;

clustrules = get_clustrules(clustrules);
 
ops0.getROIs         = getOr(db, {'getROIs'},1);
ops0.RegFile_xtime   = getOr(db, {'RegFile_xtime'},4); % number of tiffs to average when writing time averaged registered tiffs (slow)


ops = build_ops3(db(end), ops0);


clustModel     = getOr(ops, {'clustModel'}, 'standard');
neuropilSub    = getOr(ops, {'neuropilSub'}, 'surround');
splitBlocks    = getOr(ops, {'splitBlocks'}, 'none');

    
   ops1         = reg2P_kh004(ops);  % do registration
   
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

if ops.getROIs || getOr(ops, {'writeSVDroi'}, 1)
    [ops, U, Sv]    = get_svdForROI_kh(ops);
end

myfigure('Obtained EigenImage');clf;
NN=3;
for ii=1:NN^2
    mysubplot(NN,NN,ii); imagesc(U(:,:,ii));
    box off;
    axis off;
end 
 %% in case I updated python
 
 clear classes

insert(py.sys.path,int32(0),'C:\Users\hamag\.conda\envs\khtest\Lib\site-packages');
insert(py.sys.path,int32(0),'C:\home\Platex\2017Oct');
 py.sys.path
 mod = py.importlib.import_module('my_hdbscan')
py.importlib.reload(mod)

mod = py.importlib.import_module('hdbscan')
py.importlib.reload(mod)
%% 
% if you want to start from here,
% load('C:\home\ImagingExperiments\SimulatedData\SimulatedOne\20171020\1\SVDroi_SimulatedOne_20171020_plane1_ch1.mat');
% load C:\home\ImagingExperiments\B6J647\SVDroi_B6J647_CREBRC107_20170902_d2_plane1_ch1.mat;
load C:\home\ImagingExperiments\B6J647\SVDroi_B6J647_CREBRC107_20170904_d6_plane1_ch1.mat;
[Ly,Lx,nSVD]=size(U);

myfigure('Obtained EigenImage');clf;
NN=3;
for ii=1:NN^2
    mysubplot(NN,NN,ii); imagesc(U(:,:,ii));
    box off;
    axis off;
end 

npix = Ly*Lx;
U = reshape(U,[],nSVD);
% U = bsxfun(@times, U, Sv'.^.5);
U=U.*sqrt(Sv'); % > 2016b
% U=normc(U')';
%% in case to subtract neuropil (smooth components)

do_subtract_neuropil = 0; % surprizingly no effect 

if (do_subtract_neuropil)
    TileFactor = 1; % this option can be overwritten by the user
    cell_diameter = 20;
    nTiles = ceil(TileFactor * (Ly+Lx)/2 / (10 * cell_diameter)); % neuropil is modelled as nTiles by nTiles
    
    xc = linspace(1, Lx, nTiles);
    yc = linspace(1, Ly, nTiles);
    yc = yc';
    xs = 1:Lx;
    ys = 1:Ly;
    
    sigx = 4*(Lx - 1)/nTiles;
    sigy = 4*(Ly - 1)/nTiles;
    
    S = zeros(Ly, Lx, nTiles, nTiles, 'single');
    for kx = 1:nTiles
        for ky = 1:nTiles
            cosx = 1+cos(2*pi*(xs - xc(kx))/sigx);
            cosy = 1+cos(2*pi*(ys - yc(ky))/sigy);
            cosx(abs(xs-xc(kx))>sigx/2) = 0;
            cosy(abs(ys-yc(ky))>sigy/2) = 0;
            
            S(:, :,ky, kx) = cosy' * cosx;
        end
    end
    S = reshape(S, [], nTiles^2);
    S = normc(S);
    
    nBasis = nTiles^2 ;
    PixL = ones(1, Lx * Ly)';
    
    Uneu = U';
    
    Sm = bsxfun(@times, S, PixL);
    StS = Sm' * Sm;
    StU = Sm' * Uneu';
    Lam = (StS + 1e-4 * eye(nBasis)) \ StU;
    
    % recompute neuropil pixel contribution
    neuropil = Lam' * S';
    PixL = mean(bsxfun(@times, neuropil, Uneu), 1);
    PixL = bsxfun(@rdivide, PixL, mean(neuropil.^2,1)); % normalize variance to 1
    PixL = max(0, PixL);
    neuropil = bsxfun(@times, neuropil, PixL);
    U = U - neuropil'; %what's left over for cell model
    
    myfigure('Obtained EigenImage-neuropil');clf;
    NN=3;
    for ii=1:NN^2
        Utmp = reshape(U(:,ii),Ly,Lx);
        mysubplot(NN,NN,ii); imagesc(Utmp);
        box off;
        axis off;
    end
    
%     U = reshape(U,[],nSVD);
    fig_title = 'U-neuropil';
else
    fig_title = 'U-only';
end
%% construct neighbor indices
neigh = get_neighbor_index(Ly,Lx,'3x3','circular'); % I know circular boundary is not correct, but to avoid disconnected pixel, I used circular. 
%% check whether neigh vector contains neighboring index
do_check_neigh =0;

if (do_check_neigh)
    clf;
    % neigh(neigh==0)=length(neigh)+1;
    for ii=1:Ly*Lx
        A=zeros(Ly,Lx);
        ind = neigh(ii,find(neigh(ii,:)));
        A(ind)=1;
        imagesc(A); title(num2str(ii));
        pause(0.5);
    end
end
%% For each neighbor, calculate inner product of U.*(Sv.^0.5).
% To do this, construct sparse matrix which non-zero values are inner product of 
% feature vector U, but only between the 8-neighbors.

ind_i = neigh';
ind_j = repmat([1:Ly*Lx],size(neigh,2),1); % ind_j is column number.
ind_j=ind_j(:); ind_i=ind_i(:);
ind_j(ind_i==0)=[]; 
ind_i(ind_i==0)=[];

nbatch = 100000; 
nind = length(ind_i);
start_ind = 1:nbatch:nind;
end_ind   = start_ind+nbatch-1; end_ind(end)=nind;
s = zeros(1,length(ind_i));
cnt=0;
for ii=1:length(start_ind)
    index = start_ind(ii):end_ind(ii);
    s(cnt+[1:length(index)])=sum(U(ind_j(index),:).*U(ind_i(index),:),2)';
    cnt = cnt+length(index);
end

% construct sparse matrix which indicates the 8-neighbors connection.

s = max(s)-s+eps;
mins = min(s);
maxs= max(s);
s = s/(maxs-mins);

S = sparse(ind_i,ind_j,s,Ly*Lx,Ly*Lx);

%% now, how do I know I am doing right? I need to check the values of definitely correlated pixels and backgrounds
% 
% % pick the center of active cell in the image of U 
% fprintf('Select on-cell pixel: ')
% pix_on_cell =ceil(ginput(1)); 
% ind=sub2ind([Ly,Lx],pix_on_cell(2),pix_on_cell(1));
% fprintf('[%d %d] ->(%d)\n',pix_on_cell(2),pix_on_cell(1),ind);
% S(ind,:)
% 
% fprintf('Select background pixel: ')
% pix_off_cell =ceil(ginput(1)); 
% ind=sub2ind([Ly,Lx],pix_off_cell(2),pix_off_cell(1));
% fprintf('[%d %d] -> (%d)\n',pix_off_cell(2),pix_off_cell(1),ind);
% S(ind,:)

% seems OK 

%% Finally, test with hdbscan: well it did not work. Finally it worked. 
pyS = matsparse_2_pysparse(S);
 tic
%  import sys
%     sys.path.append('C:\\Users\\hamag\\.conda\\envs\\khtest\\Lib\\site-packages')
%     import hdbscan
    PA = pyargs('min_cluster_size',uint32(70),...
                'min_samples', uint32(6),...
                'gen_min_span_tree',true, ...
                'metric','precomputed',...
                'cluster_selection_method','eom',...
                'alpha',1.3);
    clusterer = py.hdbscan.HDBSCAN(PA);
    clusterer=clusterer.fit(pyS);
    toc
% clusterer=py.my_hdbscan.hdbSfit(pyS); % this is the function to do hdbscan 
 labels=nparray2mat(clusterer.labels_)+1; % python index starts from 0.
 probabilities=nparray2mat(clusterer.probabilities_);
 
 % well, it seems there is a problem in using sparse matrix especially when
 % there is a disjoint cluster. What shall I do? => set min_samples=#neighbors-1.)
%  min_samples is used to calculate the mutual_reachability.
% mutual_reachability is set to inf when #neighbors < min_samples, then it
% tends to generate disjoint data
% 

myfigure(fig_title);clf;
L=reshape(labels,Ly,Lx);
L(L==0)=length(unique(labels))+1;
% imagesc(L);
%
rand_hue  = rand(1,length(unique(labels)));
rand_hue= [rand_hue,0];
Hue = rand_hue(L);
Sat = reshape(probabilities,Ly,Lx);
Val = Sat;
IMG=hsv2rgb(cat(3,Hue,Sat,Val));
imagesc(IMG);
%% but still, some cells are merged. Why? Let's check them.
 % obtain left-top and right-bottom corner.
%  x=ceil(ginput(2));
%  J= x(1,2):x(2,2);
%  I= x(1,1):x(2,1);
%  index=reshape(1:Ly*Lx,Ly,Lx);
%  index=index(J,I);
% %  Utmp = U(index(:),:);
%  Stmp =S(index(:),index(:));
%  
%  pyS = matsparse_2_pysparse(Stmp);
%  tic
%     PA = pyargs('min_cluster_size',uint32(70),...
%                 'min_samples', uint32(2),...
%                 'gen_min_span_tree',true, ...
%                 'metric','precomputed',...
%                 'cluster_selection_method','eom',...
%                 'alpha',1);
%     clusterer2 = py.hdbscan.HDBSCAN(PA);
%     clusterer2=clusterer2.fit(pyS);
%    
%  labels2=nparray2mat(clusterer2.labels_)+1; % python index starts from 0.
%  probabilities2=nparray2mat(clusterer2.probabilities_);
%  
%  myfigure('PCA of Utmp');clf;
%  L2=reshape(labels2,length(J),length(I));
% L2(L2==0)=length(unique(labels2))+1;
% rand_hue  = rand(1,length(unique(labels2)));
% rand_hue= [rand_hue,0];
% Hue = rand_hue(L2);
% Sat = reshape(probabilities2,length(J),length(I));
% Val = Sat;
% IMG=hsv2rgb(cat(3,Hue,Sat,Val));
% imagesc(IMG);
% 
% myfigure('dendrogram of Utmp')
% slt=nparray2mat(clusterer2.single_linkage_tree_.to_numpy);
% dendrogram([slt(:,1:2)+1,slt(:,3)]); % first two columns are index starting from 0, 3rd is the distance.

%% well, dendrogram did not reveal the diffrence. Let's check the distribution of pixel.

% now, pickup the two points which you want to divide.
x = ceil(ginput(2));

%
% if iclust1~=iclust2
%     error('Please select the same ROI to split');
% end
index=reshape(1:Ly*Lx,Ly,Lx);
indY=[];
indX=[];
indAll = [];
seed_ind = zeros(size(x,1),1);
seed_ind_in_indAll = zeros(size(x,1),1);
iclust=zeros(size(x,1),1);
Val = zeros(Ly,Lx);
for ii=1:size(x,1);
    tmp=L(x(ii,2),x(ii,1));
    fprintf('iclust%d=%d\n ',ii,tmp);
    if any(iclust==tmp)
        % detected same cluster
        fprintf('Detected the same cluster, skip!\n');
        iclust(ii)=[];
        continue;
    else
        iclust(ii) = tmp;
    end
    ind=find(L==iclust(ii));
    % 
    
    % get the square region 
    [index_y,index_x]=ind2sub([Ly,Lx],ind);
    indY=cat(1,indY,index_y);
    indX=cat(1,indX,index_x);
    indAll{ii} = ind;
    Val = Val | L==iclust(ii);
    
    seed_ind(ii)=index(x(ii,2),x(ii,1));
    seed_ind_in_indAll(ii)=find(ind==seed_ind(ii));
end
J=min(indY):max(indY);
I=min(indX):max(indX);
myfigure('Utmp plot');clf;
subplot(2,2,1);

Hue = rand_hue(L);
Sat = reshape(probabilities,Ly,Lx);
IMG=hsv2rgb(cat(3,Hue,Sat,Val));
imagesc(IMG(J,I,:));


%%
all_index = cat(1,indAll{:});
Utmp = U(all_index,:);
[coef,score] = pca(Utmp);
subplot(2,2,2);cla;
plotH=[];
ind=0;
Col = {'r','g','b'}; 
for ii=1:length(indAll)
    ind=ind+[1:length(indAll{ii})];
    plotH(ii)=plot3(score(ind,1),score(ind,2),score(ind,3),[Col{ii},'.']);hold on;
    indtmp=ind(seed_ind_in_indAll(ii));
    [score(indtmp,1), score(indtmp,2), score(indtmp,3)]
    plot3(score(indtmp,1), score(indtmp,2), score(indtmp,3),[Col{ii},'o']);
        ind = ind(end); 
%       pause
end

grid on;
%% next, plot the result of k-means

Nlabel = 2;
 idx = kmeans(Utmp,Nlabel);
 subplot(2,2,3);cla;
 for ii=1:Nlabel
     ind=find(idx==ii);
       plotH(ii)=plot3(score(ind,1),score(ind,2),score(ind,3),[Col{ii},'.']);hold on;
 end
grid on;

%% recover the divided image 

subplot(2,2,4);cla;

Hue = rand_hue(L);
Sat = reshape(probabilities,Ly,Lx);

Col = [0.9,0.4,0.8];
% Val = L~=0;

for ii=1:Nlabel
    original_ind = all_index(idx==ii);
    Hue(original_ind)=Col(ii);
%     Sat(original_ind)=1;
end

IMG=hsv2rgb(cat(3,Hue,Sat,Val));
image(IMG(J,I,:));
%% let's see how dbscan can divide this 
Spind = L(J,I); Spind = Spind(:);

%% now compare different parameters

alpha = [1,1.3, 1.6];
min_cluster_size = [50,100,150];
min_samples = [2,4,6];



for ms = min_samples
    
    myfigure(sprintf('min_samples=%d',ms));clf;
    cnt=0;
    for mc  = min_cluster_size
        for aa=alpha
            cnt=cnt+1;
            subplot(length(alpha),length(min_cluster_size),cnt);
            titletxt=sprintf('alpha=%3.2f, min_cluser_size=%d, min_samples=%d\n',aa,mc,ms);
            PA = pyargs('min_cluster_size',uint32(mc),...
                'min_samples', uint32(ms),...
                'gen_min_span_tree',true, ...
                'metric','precomputed',...
                'cluster_selection_method','leaf',...
                'alpha',aa);
            clusterer = py.hdbscan.HDBSCAN(PA);
            clusterer=clusterer.fit(pyS);
            labels=nparray2mat(clusterer.labels_)+1; % python index starts from 0.
            probabilities=nparray2mat(clusterer.probabilities_);
            L=reshape(labels,Ly,Lx);
            imagesc(L);
            title(titletxt);
            drawnow;
        end
    end
end

   

%  Stril = tril(S);
% [I,J,Stril]=find(Stril);
%% different approach 2
pyS = matsparse_2_pysparse(S);
% py.sys.path.append('C:\\Users\\hamag\\.conda\\envs\\khtest\\Lib\\site-packages')

  lil_matrix = pyS.tolil();
  mutual_reachability_ = py.hdbscan.hdbscan_.sparse_mutual_reachability(lil_matrix, pyargs('min_points',2));
%   n=py.scipy.sparse.csgraph.connected_components(mutual_reachability_,pyargs('directed',false,'return_labels',true));
%   if n>1 error('min_points are too small, some clusters are completely separated'); end
  sparse_min_spanning_tree = py.scipy.sparse.csgraph.minimum_spanning_tree(mutual_reachability_);
  
   nonzeros = sparse_min_spanning_tree.nonzero()
    nonzero_vals = sparse_min_spanning_tree[nonzeros]
    min_spanning_tree = np.vstack(nonzeros + (nonzero_vals,)).T
    
%% and the following is the standard suite2P code for clustering.
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
