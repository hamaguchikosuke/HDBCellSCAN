function Fcell = get_signals_and_neuropil_kh(ops, PlaneChString)

if ~isfield(ops,'zoomMicro')
    ops.zoomMicro=1; %fixed zoomMicro
end
if ~isfield(ops,'inNeurop')
    ops.inNeurop=3; %fixed inner diameter of the neuropil mask donut
end
if ~isfield(ops,'outNeurop')
    ops.outNeurop=45; %radius of Neuropil fixed at 45um
end
if ~isfield(ops,'microID')
%     ops.microID='b'; %microscope identity. 'b': b-scope. 'm': 
    ops.microID='o'; % olympus
end
if ~isfield(ops, 'processed_date') % use processed data in F_...._proc (generated in gui2P)
    UseProcFile = 0;
else
    UseProcFile = 1;
end

if ~isfield(ops, 'useSVD') % redo calculation of signal and neuropil based on
                         % saved SVD components instead of temporary
                         % bin-file
    ops.useSVD = 0;
end
if ~isfield(ops, 'getSignal') % calculate and save neural signal
    ops.getSignal = 1;
end
if ~isfield(ops, 'getNeuropil') % calculate and save neuropil
    ops.getNeuropil = 1;
end
if ~isfield(ops, 'newFile') % save new file '<name>_new.mat', otherwise
                            % existing file is overwritten
    ops.newFile = 0;
end
%%
if UseProcFile
    LoadName = ops.ProcFileName;
    SaveName = LoadName;
else
    LoadName = ops.ROISaveFile;
    SaveName = strrep(LoadName,'ROI','Fsig');
end
if exist(LoadName,'file')
%     fprintf('Loading %s\n',RoiFile);
else
    errordlg(sprintf('%s does not exist!',LoadName));
end

[SavePath,SaveName,Ext]=fileparts(SaveName);

if ops.newFile
    SaveName = [SaveName,'_new'];
end
data = load(LoadName);
% if UseProcFile == 1
%     data = data.dat;
% end
%%
ops = update_ops(ops,data.ops);
% ops = data.ops;
Nk       = numel(data.stat); % all ROIs including parents
% if UseProcFile == 1
%     Nk_parents = find(isinf(data.cl.Mrs), 1, 'last');
% else
%     Nk_parents = 0; % no parents in HDBCellScan
% end

if UseProcFile == 1
    useCells = find(data.cl.selected);
else
    useCells = find([data.stat.igood]);
end

[LyU, LxU] = size(ops.mimg);
LyR=length(ops.yrange);
LxR=length(ops.xrange);

% allField=zeros(LyR,LxR);

cellFields=zeros(length(useCells),LyR,LxR);
% figure(1001);clf;
for ii=1:length(useCells)
    ipix=data.stat(useCells(ii)).ipix;
%     temp=zeros(LyR,LxR);
%     temp(ipix)=1;
%     imagesc(temp); title(num2str(ii)); pause; 
    cellFields(ii,ipix)=1; % fast data filling 
end
 allField=squeeze(sum(cellFields,1));
ops.totPixels=LxU;

um2pix=infoPixUm(ops.totPixels,ops.zoomMicro,ops.microID);
xPU=um2pix.xPU;
yPU=um2pix.yPU;

%%% generate neuropil masks %%%
if ops.getNeuropil
    neuropMasks=createNeuropilMasks(cellFields,allField,xPU,yPU,ops);

    mCell=0;
    for k=useCells(:)'
        mCell=mCell+1;
        tmp=squeeze(neuropMasks(mCell,:,:));
        data.stat(k).ipix_neuropil=find(tmp);
    end
end

%% get signals and neuropil
if ops.useSVD == 0
    nimgbatch = 2000;
    ix = 0;
    fclose all;
    
%     fid = fopen(ops.RegFile, 'r');
    
    [SearchTiffString,TiffPaths]=get_tiff_location_from_ops_20171102(ops,PlaneChString);
    TiffFiles = {};
    for ii=1:length(TiffPaths)
        TiffFiletmp=dir(fullfile(TiffPaths{ii},SearchTiffString));
        TiffFiles = cat(2,TiffFiles,fullfile(TiffPaths{ii},{TiffFiletmp.name}));
    end
    % multipage tiff class
    MTiff = mlttiff(TiffFiles);
    
    tic
    F    = NaN(Nk, MTiff.TotalNSeries, 'single');
    Fneu = NaN(Nk, MTiff.TotalNSeries, 'single');
    
    StartInd=MTiff.CumNSeries(1:end-1)+1;
    EndInd = MTiff.CumNSeries(2:end);

    for ii=1:length(MTiff.nSeries)
%         mov = fread(fid,  LyU*LxU*nimgbatch, '*int16');
        index = StartInd(ii):EndInd(ii);
        mov = MTiff.get_stacks(index);
      
%         mov = reshape(mov, LyU, LxU, []);
        mov = mov(ops.yrange, ops.xrange, :);
        mov = single(mov);
        NT= size(mov,3);
        
        mov = single(reshape(mov, [], NT));
        
        for k = 1:Nk
            if ops.getSignal
                ipix = data.stat(k).ipix;
                if ~isempty(ipix)
                    % F(k,ix + (1:NT)) = stat(k).lambda' * data(ipix,:);
                    F(k,index) = mean(mov(ipix,:), 1);
                end
            end
            if ops.getNeuropil
                ipix_neuropil= data.stat(k).ipix_neuropil;
                if ~isempty(ipix_neuropil)
                    Fneu(k,index) = mean(mov(ipix_neuropil,:), 1);
                end
            end
        end
        
        
        fprintf('Frame %d done in time %2.2f \n', index(end), toc)
        
    end
    % F = F(:, 1:ix);
    if cumsum(ops.Nframes)~=MTiff.TotalNSeries
        QText=sprintf('The number of frames in image registration is %d, but number of frames found in Tiff file was %d. Can I limit the total frame to %d? ',...
            cumsum(ops.Nframes),MTiff.TotalNSeries,MTiff.TotalNSeries);
        button = questdlg(QText,'Frame number mismatch');
        switch button
            case 'Yes'
%                 NframeOps = MTiff.TotalNSeries;
                ops.Nframes = MTiff.TotalNSeries;
            case {'No', 'Cancel'}
                disp('User cancelled.');
                return;
        end
    end

    csumNframes = [0 cumsum(ops.Nframes)];
    Fcell = cell(1, length(ops.Nframes));
    FcellNeu = cell(1, length(ops.Nframes));
    for ii = 1:length(ops.Nframes)
        Fcell{ii} 	= F(:, csumNframes(ii) + (1:ops.Nframes(ii)));
        FcellNeu{ii} = Fneu(:, csumNframes(ii) + (1:ops.Nframes(ii)));
    end
else % ops.useSVD == 1
    ind = strfind(filenames{1}, '_Nk');
    fname = ['SVD' filenames{1}(2:ind-1) '.mat'];
    svd = load(fullfile(ops.ResultsSavePath, fname), 'U', 'Vcell');
    U1 = reshape(svd.U, [], size(svd.U,3));
    Ucell = zeros(length(useCells), size(U1,2));
    UcellNeu = zeros(length(useCells), size(U1,2));
    Fcell = cell(size(svd.Vcell));
    FcellNeu = cell(size(svd.Vcell));
    for iCell = 1:length(useCells)
        if ops.getSignal
            ipix = data.stat(useCells(iCell)).ipix;
            if ~isempty(ipix)
                Ucell(iCell,:) = mean(U1(ipix, :), 1);
            end
        end
        if ops.getNeuropil
            ipix_neuropil= data.stat(useCells(iCell)).ipix_neuropil;
            if ~isempty(ipix_neuropil)
                UcellNeu(iCell,:) = mean(U1(ipix_neuropil,:), 1);
            end
        end
    end
    for iExp = 1:length(svd.Vcell)
        if ops.getSignal
            F = NaN(Nk, size(svd.Vcell{iExp}, 2), 'single');
            F(useCells,:) = Ucell * svd.Vcell{iExp};
            Fcell{iExp} = F;
        end
        if ops.getNeuropil
            F = NaN(Nk, size(svd.Vcell{iExp}, 2), 'single');
            F(useCells,:) = UcellNeu * svd.Vcell{iExp};
            FcellNeu{iExp} = F;
        end
    end
end

if UseProcFile
    if ops.getSignal
        data.F.Fcell = Fcell;
    end
    if ops.getNeuropil
        data.F.FcellNeu = FcellNeu;
    end
else
    if ops.getSignal
        data.Fcell = Fcell;
    end
    if ops.getNeuropil
        data.FcellNeu = FcellNeu;
    end
end

dat = data;
dat.opsNpil = ops;
FullSaveName = fullfile(SavePath,SaveName);
fprintf('Saving signal in %s\n',FullSaveName);
save(FullSaveName, '-struct', 'dat');


