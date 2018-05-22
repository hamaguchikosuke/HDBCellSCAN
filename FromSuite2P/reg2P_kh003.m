function ops1 = reg2P_kh003(ops)
%% Align GCaMP and RCaMP data all together.
% nchannels:
%
numPlanes = length(ops.planesToProcess);

nplanes             = getOr(ops, {'nplanes'}, 1);
nchannels           = getOr(ops, {'nchannels'}, 1);

Ca_channel          = getOr(ops, {'Ca_channel'}, 1);
Ref_channel         = getOr(ops, {'Align_channel'}, 2);

if length(Ca_channel)>nchannels
    error('Length of Ca_channel must be larger than nchannels')
end

BiDiPhase            = getOr(ops, {'BiDiPhase'}, 0);
LoadRegMean         = getOr(ops, {'LoadRegMean'}, 0);
ops.RegFileBinLocation = getOr(ops, {'RegFileBinLocation'}, []);
ops.splitFOV           = getOr(ops, {'splitFOV'}, [1 1]);


ops.smooth_time_space = getOr(ops, 'smooth_time_space', []);
fs = ops.fsroot;

%% find the mean frame after aligning a random subset
ntifs = sum(cellfun(@(x) numel(x), fs));
nfmax = ceil(ops.NimgFirstRegistration/ntifs);
if nfmax>=2000
    nfmax = 1999;
end

Info0 = imfinfo(fs{1}(1).name);
Ly = Info0(1).Height;
Lx = Info0(1).Width;
ops.Lx = Lx;
ops.Ly = Ly;

[xFOVs, yFOVs] = get_xyFOVs(ops);
ops.NumSplitViews=prod(ops.splitFOV);
ops.xFOVs = xFOVs;
ops.yFOVs = yFOVs;

if ops.doRegistration
    
    IMG = get_sparse_subset_images(ops.SubDirs,fs,nchannels,nplanes,Ref_channel,ops);
    % IMG : Height x Width x nPlane x nSeries
    
    ops1 = cell(numPlanes, ops.NumSplitViews,nchannels);
    % First, process the reference image.
    for ii = 1:numPlanes
        for jj = 1:ops.NumSplitViews
            ops1{ii,jj,Ref_channel} = align_iterative(single(squeeze(IMG(yFOVs(:,jj),xFOVs(:,jj),ops.planesToProcess(ii),:))), ...
                ops);
            % use ops.dsprealign to align other images
            ops1{ii,jj,Ref_channel}.RefImg=ops1{ii,jj,Ref_channel}.mimg;
            ops1{ii,jj,Ref_channel}.IsRefChannel= 1;
            ops1{ii,jj,Ref_channel}.PlaneID = ii;
            ops1{ii,jj,Ref_channel}.ViewID  = jj;
            ops1{ii,jj,Ref_channel}.ChannelID  = Ref_channel;
            
        end
    end
    
    
    % In case multi-color Ca imaging.
    OtherCaChannels = setdiff(Ca_channel ,Ref_channel);
    
    for cc=1:length(OtherCaChannels)
        IMG = get_sparse_subset_images(ops.SubDirs,fs,nchannels,nplanes,OtherCaChannels(cc),ops);
        % IMG : Height x Width x nPlane x nSeries
        
        
        for ii = 1:numPlanes
            for jj = 1:ops.NumSplitViews
                ops1{ii,jj,OtherCaChannels(cc)} = ops1{ii,jj,Ref_channel}; % just copy contents from Ref_channel. 
                ops1{ii,jj,OtherCaChannels(cc)}.RefImg = ...
                    mean(register_movie(single(squeeze(IMG(yFOVs(:,jj),xFOVs(:,jj),ops.planesToProcess(ii),:))), ...
                    ops, ops1{ii,jj,Ref_channel}.dsprealign),3);
                ops1{ii,jj,OtherCaChannels(cc)}.IsRefChannel= 0;
                ops1{ii,jj,OtherCaChannels(cc)}.PlaneID = ii;
                ops1{ii,jj,OtherCaChannels(cc)}.ViewID  = jj;
                ops1{ii,jj,OtherCaChannels(cc)}.ChannelID  = OtherCaChannels(cc);
            
                % use ops.dsprealign to align other images
            end
        end
    end
    
    if ops.showTargetRegistration
        
        for cc=1:length(Ca_channel)
            h=myfigure(sprintf('TargetImage:Ch%d',Ca_channel(cc)));
            set(h,'position', [900 50 900 300*numPlanes])
            ax = ceil(sqrt(numel(ops1)/2));
            i0 = 0;
            for ii = 1:numPlanes
                for jj = 1:ops.NumSplitViews
                    i0 = i0+1;
                    subplot(ax,2*ax,i0)
                    imagesc(ops1{ii,jj,cc}.RefImg);
                    colormap('gray')
                    title(sprintf('Registration for Ch=%d, plane %d, mouse %s, date %s', ...
                        Ca_channel(cc), ii, ops.mouse_name, ops.date))
                end
            end
        end
        
        drawnow
    end
    
    clear IMG
else
    for ii = 1:numPlanes
        ops1{ii} = ops;
        ops1{ii}.mimg = zeros(Ly, Lx);
        ops1{ii}.Ly   = Ly;
        ops1{ii}.Lx   = Lx;
    end
end
%%
fid = cell(numPlanes, ops.NumSplitViews,length(Ca_channel));
for cc=1:length(Ca_channel)
    for ii = 1:numPlanes
        for jj = 1:ops.NumSplitViews
            ops1{ii,jj,cc}.RegFile = fullfile(ops.RegFileRoot, sprintf('tempreg_plane%d_ch%d.bin', ii + (jj-1)*numPlanes,cc));
            regdir = fileparts(ops1{ii,jj,cc}.RegFile);
            if ~exist(regdir, 'dir')
                mkdir(regdir);
            end
            
            % open bin file for writing
            fid{ii,jj,cc}              = fopen(ops1{ii,jj,cc}.RegFile, 'w');
            ops1{ii,jj,cc}.DS          = [];
            ops1{ii,jj,cc}.CorrFrame   = [];
            ops1{ii,jj,cc}.mimg1       = zeros(ops1{ii}.Ly, ops1{ii}.Lx);
            
        end
    end
end


%% Then by using the reference frame, align all the frames.
tic

for kk = 1:length(fs) % For each SubDir
    for ii = 1:numel(ops1)
        ops1{ii}.Nframes(kk)  = 0;
    end
    
    iplane0 = 1:1:ops.nplanes;
    LastnFr = 0;
    for ff = 1:length(fs{kk}) % For each Tiff file
%         iplane0 = mod(iplane0-1, numPlanes) + 1;
        iplane0 = circshift(iplane0,-mod(LastnFr/nchannels,ops.nplanes));
        fprintf('Copying %s to local (%s)\n',fs{kk}(ff).name,ops.temp_tiff);
        copyfile(fs{kk}(ff).name,ops.temp_tiff);
        
        nFr = nFrames(ops.temp_tiff);
        if mod(nFr, nchannels) ~= 0
            fprintf('  WARNING: number of frames in tiff (%d) is NOT a multiple of number of channels!\n', ff);
        end
        nChunk = 6000*nchannels;
        nSteps = ceil(nFr/nChunk);
        StartFrames =1+[0:nChunk:nChunk*(nSteps-1)];
        EndFrames   = [nChunk:nChunk:nChunk*nSteps];
        EndFrames(end)=min(nFr,EndFrames(end));
        
        for ss=1:nSteps
            % Because the data is read as a chunk, starting position is not
            % always plane 1. 
            
            for cc=[Ref_channel,OtherCaChannels] % start from Ref_channel, after making ops1.dsall
                if ops.doRegistration
                 
                    scan_index = (StartFrames(ss)+cc-1):nchannels:EndFrames(ss);
  
                    % Calculate Reference Image for each channel, based on the movement in reference channel.
                    
                    % Reference Image: mean registered image of its own channel
                    
                    % reference channel: channel within in which movements are used as
                    % a reference for other channels. 
                    
                    if cc==Ref_channel
                        [ops1,data,dsall] = partial_reg(ops,ops1,cc,scan_index,iplane0);
                    else        
                        data = BigTiffReader(ops.temp_tiff, scan_index);
                    end
                    
                    
                    
                    % ix0 = StartFrames(ss)-1;
                    ix0 = 0;
                    Nbatch = 1000;
                    dreg = zeros(size(data), class(data));
                    
                    while ix0<size(data,3)
                        indxr = ix0 + (1:Nbatch);
                        indxr(indxr>size(data,3)) = [];
%                         ds_index = indxr+StartFrames(ss)-1;
                        for ll = 1:ops.NumSplitViews
%                             dreg(yFOVs(:,ll), xFOVs(:,ll), indxr)        = ...
%                                 register_movie(data(yFOVs(:,ll), xFOVs(:,ll), indxr), ...
%                                 ops1{1,ll,cc}, ops1{1,ll,Ref_channel}.DS(ds_index,:,ll));
                               dreg(yFOVs(:,ll), xFOVs(:,ll), indxr)        = ...
                                register_movie(data(yFOVs(:,ll), xFOVs(:,ll), indxr), ...
                                ops1{1,ll,cc}, dsall(indxr,:,ll));
                        end
                        ix0 = ix0 + Nbatch;
                    end
                    
                else   %----------  NO registration -------------%
                    scan_index = (StartFrames(ss)+cc-1):nchannels:EndFrames(ss);
                    dreg = BigTiffReader(ops.temp_tiff, scan_index);  % dreg is a chunk of registered image.
                end
                
                % write dreg to bin file+
                for ii = 1:numPlanes
                    ifr0 = iplane0(ops.planesToProcess(ii));
                    indframes = ifr0:nplanes:size(data,3);
                    for ll = 1:ops.NumSplitViews
                        dwrite = dreg(yFOVs(:,ll),xFOVs(:,ll),indframes);
                        datatype = class(data);
                        fwrite(fid{ii,ll,cc}, dwrite, datatype);
                        
                        ops1{ii,ll,cc}.Nframes(kk) = ops1{ii,ll,cc}.Nframes(kk) + size(dwrite,3);
                        ops1{ii,ll,cc}.mimg1 = ops1{ii,ll,cc}.mimg1 + sum(dwrite,3);
                    end
                end
                
                fprintf('Chunk %d (%d/%d) \n', ss,EndFrames(ss),nFr)
                
            end
            if rem(ff,5)==1
                fprintf('Set %d, tiff %d done in time %2.2f \n', kk, ff, toc)
            end
            LastnFr = nFr;
            
%             iplane0 = iplane0 - nFr/nchannels;
        end
        
    end
end 
   %%
    for ii = 1:numel(ops1)
        ops1{ii}.mimg1 = ops1{ii}.mimg1/sum(ops1{ii}.Nframes);
    end
    
    %% register all the Ca images. 
    for ii = 1:numel(ops1)
        fclose(fid{ii});
        if ~any(ops1{ii}.ChannelID == Ca_channel)
            continue; % this channel is not Ca-image. Don't bother to save the registered data. 
        else
            
            fprintf('Loading %s...\n',ops1{ii}.RegFile);
            fid{ii}           = fopen(ops1{ii}.RegFile, 'r');
            
            if ~isempty(ops.RegFileTiffLocation)
                ops1{ii} = write_reg_to_tiff2(fid{ii}, ops1{ii}, ii,datatype);
                write_reg_tave_to_tiff2(fid{ii}, ops1{ii}, ii,ops.RegFile_xtime,datatype);
            end
            
            if ~isempty(ops.nimgbegend) && ops.nimgbegend>0
                ops1{ii} = getBlockBegEnd(fid{ii}, ops1{ii}); % get mean of first and last frames in block (to check for drift)
            end
            
            if ~isempty(ops.RegFileBinLocation)
                folder = fullfile(ops1{ii}.RegFileBinLocation, ops1{ii}.mouse_name, ...
                    ops1{ii}.date);
                if ~exist(folder, 'dir')
                    mkdir(folder)
                end
                [~,regfile,~]=fileparts(ops1{ii}.RegFile);
                regfile = strrep(regfile,'tempreg_','');
                fidCopy = fopen(fullfile(folder, sprintf('%s.bin', regfile)), 'w');
                sz = ops1{ii}.Lx * ops1{ii}.Ly;
                parts = ceil(sum(ops1{ii}.Nframes) / 2000);
                for p = 1:parts
                    toRead = 2000;
                    if p == parts
                        toRead = sum(ops1{ii}.Nframes) - 2000 * (parts-1);
                    end
                    data = fread(fid{ii},  sz*toRead, '*int16');
                    fwrite(fidCopy, data, class(data));
                end
                fclose(fidCopy);
            end
            fclose(fid{ii});
        end
    end
    %% before saving ops, distributed the mean image data in each ops{}.
    
     for ii = 1:numel(ops1)
         
     end
        
    %%
    % compute xrange, yrange
    for ii = 1:numel(ops1)
        if ops.doRegistration && (ops1{ii}.IsRefChannel)
            minDs = min(ops1{ii}.DS(2:end, [1 2]), [], 1);
            maxDs = max(ops1{ii}.DS(2:end, [1 2]), [], 1);
            disp([minDs(1) maxDs(1) minDs(2) maxDs(2)]);
            if BiDiPhase>0
                maxDs(2) = max(1+BiDiPhase, maxDs(2));
            elseif BiDiPhase<0
                minDs(2) = min(BiDiPhase, minDs(2));
            end
            
            ops1{ii}.yrange = ceil(maxDs(1)):floor(ops1{ii}.Ly+minDs(1));
            ops1{ii}.xrange = ceil(maxDs(2)):floor(ops1{ii}.Lx+minDs(2));
        else
            ops1{ii}.yrange = 1:Ly;
            ops1{ii}.xrange = 1:Lx;
        end
        
        savepath = sprintf('%s/', ops.ResultsSavePath);
        
        if ~exist(savepath, 'dir')
            mkdir(savepath)
        end
          [~,regfile,~]=fileparts(ops1{ii}.RegFile);
            regfile = strrep(regfile,'tempreg_','');
            
        ops = ops1{ii};
        save(sprintf('%s/regops_%s_%s_%s.mat', ops.ResultsSavePath, ...
            ops.mouse_name, ops.date, regfile),  'ops')
    end

     
%     save(sprintf('%s/F_%s_%s_All.mat', ops.ResultsSavePath, ...
%       ops.mouse_name, ops.date, iplane), 'ops')

    %%
    function [ops1,data,dsall] = partial_reg(ops,ops1,ch_to_align,scan_index,iplane0)
    
    NumPlanes = length(ops.planesToProcess);
    nplanes             = getOr(ops, {'nplanes'}, 1);

    data = BigTiffReader(ops.temp_tiff, scan_index);
    
    BiDiPhase            = getOr(ops, {'BiDiPhase'}, 0);
    
    if BiDiPhase
        data = ResonanceImagingPhaseShifter(data,BiDiPhase);
    end
    
    % get the registration offsets
    dsall = zeros(size(data,3), 2, ops.NumSplitViews);
    for ii = 1:NumPlanes
        ifr0 = iplane0(ops.planesToProcess(ii));
        indframes = ifr0:nplanes:size(data,3);
        
        for ll = 1:ops.NumSplitViews
            dat = data(ops.yFOVs(:,ll),ops.xFOVs(:,ll),indframes);
            if ~isempty(ops.smooth_time_space)
                dat = smooth_movie(dat, ops);
            end
            [ds, Corr]  = registration_offsets(dat, ops1{ii,ll,ch_to_align}, 0,ops1{ii,ll,ch_to_align}.RefImg);
            dsall(indframes,:, ll)  = ds;
%             % collect ds; wait, I think it is not necessary. 
%             if ff==1
%                 ds(1,:,:) = 0;
%             end
            ops1{ii,ll,ch_to_align}.DS          = cat(1, ops1{ii,ll,ch_to_align}.DS, ds);
            ops1{ii,ll,ch_to_align}.CorrFrame   = cat(1, ops1{ii,ll,ch_to_align}.CorrFrame, Corr);
        end
    end
    
    % added by KH 20161122 to suppress extremely large movement.
    MovementShift = sqrt(sum(ops1{ii,ll,ch_to_align}.DS.^2,2));
    TooBigMovementIndex = find(MovementShift>ops.MaxMovementPixel);
    for ii=1:TooBigMovementIndex
        if TooBigMovementIndex(ii)>1 || TooBigMovementIndex(ii)<length(ops1{ii,ll,ch_to_align}.DS)
            fprintf('!');
            ops1{ii,ll,ch_to_align}.DS(TooBigMovementIndex(ii))= 0.5*( ops1{ii,ll,cc,ch_to_align}.DS(TooBigMovementIndex(ii-1)) +  ops1{ii,ll,ch_to_align}.DS(TooBigMovementIndex(ii+1)))
        end
    end
