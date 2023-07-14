function ops1 = reg2P_kh007(ops)
%% Align GCaMP and RCaMP data all together.
% 006 updates: GrinLens mode is added.
% 
% by Kosuke Hamaguchi
numPlanes = length(ops.planesToProcess);

nplanes             = getOr(ops, {'nplanes'}, 1);
nchannels           = getOr(ops, {'nchannels'}, 1);

Ca_channel          = getOr(ops, {'Ca_channel'}, 1);
Ref_channel         = getOr(ops, {'Align_channel'}, 2);

ops.RegShape            = getOr(ops, {'RegShape'},'Square');
fprintf('RegShape=%s\n',ops.RegShape)

if length(Ca_channel)>nchannels
    error('Length of Ca_channel must be larger than nchannels')
end

dobidi             = getOr(ops, {'dobidi'}, 1); % compute bidiphase?
% if set to a value by user, do not recompute
if isfield(ops, 'BiDiPhase')
    dobidi         = 0;
end 

BiDiPhase            = getOr(ops, {'BiDiPhase'}, 0);
LoadRegMean         = getOr(ops, {'LoadRegMean'}, 0);
ops.RegFileBinLocation = getOr(ops, {'RegFileBinLocation'}, []);
ops.splitFOV           = getOr(ops, {'splitFOV'}, [1 1]);


ops.smooth_time_space = getOr(ops, 'smooth_time_space', []);
fs = ops.fsroot;

if isempty(fs{1})
    error('No .tif files in %s\n',fullfile(ops.RootDir,ops.SubDirs{1}));
end

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

ops.RegRange_Xind= getOr(ops, {'RegRange_Xind'},1:ops.Lx);
ops.RegRange_Yind= getOr(ops, {'RegRange_Yind'},1:ops.Ly);

switch ops.RegShape
    case 'Square'
        [xFOVs, yFOVs] = get_xyFOVs(ops);
        ops.NumSplitViews=prod(ops.splitFOV);
        ops.xFOVs = xFOVs;
        ops.yFOVs = yFOVs;

    case {'GrinLens','TopLeft','TopCenter','TopRight','CenterLeft'}
        [xFOVs, yFOVs] = get_xyFOVs(ops,ops.RegShape);
        ops.NumSplitViews=1;
        ops.xFOVs = xFOVs;
        ops.yFOVs = yFOVs;

    otherwise
        error();
end
if ops.doRegistration
    
    IMG = get_sparse_subset_images(ops.SubDirs,fs,nchannels,nplanes,Ref_channel,ops);
     if (ops.IMGLogScaling), IMG = log(single(IMG)+1); end
    % IMG : Height x Width x nPlane x nSeries
    if dobidi
        ops.BiDiPhase = BiDiPhaseOffsets(IMG);
        fprintf('Estimated Phase difference=%d\n',ops.BiDiPhase)
    end
    ops1 = cell(numPlanes, ops.NumSplitViews,nchannels);
    % First, process the reference image.
    for ii = 1:numPlanes
        for jj = 1:ops.NumSplitViews
            ops1{ii,jj,Ref_channel} = align_iterative_for_HDBCellSCAN(single(squeeze(IMG(yFOVs(:,jj),xFOVs(:,jj),ops.planesToProcess(ii),:))), ...
                ops);
            % use ops.dsprealign to align other images
            ops1{ii,jj,Ref_channel}.RefImg=ops1{ii,jj,Ref_channel}.mimg;
            ops1{ii,jj,Ref_channel}.IsRefChannel= 1;
            ops1{ii,jj,Ref_channel}.PlaneID = ii;
            ops1{ii,jj,Ref_channel}.ViewID  = jj;
            ops1{ii,jj,Ref_channel}.ChannelID  = Ref_channel;
            
        end
    end
    
    
    % In case multi-channel imaging.
%     OtherCaChannels = setdiff(Ca_channel ,Ref_channel);
     RemainingChannels = setdiff([1:ops.nchannels] ,Ref_channel);
    
    for cc=1:length(RemainingChannels)
        IMG = get_sparse_subset_images(ops.SubDirs,fs,nchannels,nplanes,RemainingChannels(cc),ops);
        % IMG : Height x Width x nPlane x nSeries
        %<ToDO>: apply log scaling 
        if (ops.IMGLogScaling), IMG = log(single(IMG)+1); end
        for ii = 1:numPlanes
            for jj = 1:ops.NumSplitViews
                ops1{ii,jj,RemainingChannels(cc)} = ops1{ii,jj,Ref_channel}; % just copy contents from Ref_channel. 
                ops1{ii,jj,RemainingChannels(cc)}.RefImg = ...
                    mean(register_movie(single(squeeze(IMG(yFOVs(:,jj),xFOVs(:,jj),ops.planesToProcess(ii),:))), ...
                    ops, ops1{ii,jj,Ref_channel}.dsprealign),3);
                ops1{ii,jj,RemainingChannels(cc)}.IsRefChannel= 0;
                ops1{ii,jj,RemainingChannels(cc)}.PlaneID = ii; 
                ops1{ii,jj,RemainingChannels(cc)}.ViewID  = jj;
                ops1{ii,jj,RemainingChannels(cc)}.ChannelID  = RemainingChannels(cc);
            
                % use ops.dsprealign to align other images
            end
        end
    end
    
    if ops.showTargetRegistration
        
        for cc=1:ops.nchannels
            h=myfigure(sprintf('TargetImage:Ch%d',cc));
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
                        cc, ii, ops.mouse_name, ops.date))
                end
            end
        end
        
        drawnow
    end
    
    clear IMG
else % no registration. Just add an average of a sparse subset of images as mimg 
    IMG = get_sparse_subset_images(ops.SubDirs,fs,nchannels,nplanes,Ref_channel,ops);
    % IMG : Height x Width x nPlane x nSeries
    
    ops1 = cell(numPlanes, ops.NumSplitViews,nchannels);
    % First, process the reference image.
    for ii = 1:numPlanes
        for jj = 1:ops.NumSplitViews
            ops1{ii,jj,Ref_channel} = ops;
            ops1{ii,jj,Ref_channel}.mimg = mean(single(squeeze(IMG(yFOVs(:,jj),xFOVs(:,jj),ops.planesToProcess(ii),:))),3);
            % use ops.dsprealign to align other images
            ops1{ii,jj,Ref_channel}.RefImg=ops1{ii,jj,Ref_channel}.mimg;
            ops1{ii,jj,Ref_channel}.IsRefChannel= 1;
            ops1{ii,jj,Ref_channel}.PlaneID = ii;
            ops1{ii,jj,Ref_channel}.ViewID  = jj;
            ops1{ii,jj,Ref_channel}.ChannelID  = Ref_channel;
           
            ops1{ii,jj,Ref_channel}.Ly   = Ly;
            ops1{ii,jj,Ref_channel}.Lx   = Lx;
      
        end
    end
    
    
    % In case multi-color Ca imaging.
    RemainingChannels = setdiff(Ca_channel ,Ref_channel);
    
    for cc=1:length(RemainingChannels)
        IMG = get_sparse_subset_images(ops.SubDirs,fs,nchannels,nplanes,RemainingChannels(cc),ops);
        % IMG : Height x Width x nPlane x nSeries
        
        
        for ii = 1:numPlanes
            for jj = 1:ops.NumSplitViews
                ops1{ii,jj,RemainingChannels(cc)} = ops1{ii,jj,Ref_channel}; % just copy contents from Ref_channel.
                ops1{ii,jj,RemainingChannels(cc)}.RefImg = ...
                    mean(single(squeeze(IMG(yFOVs(:,jj),xFOVs(:,jj),ops.planesToProcess(ii),:))),3);
                ops1{ii,jj,RemainingChannels(cc)}.IsRefChannel= 0;
                ops1{ii,jj,RemainingChannels(cc)}.PlaneID = ii;
                ops1{ii,jj,RemainingChannels(cc)}.ViewID  = jj;
                ops1{ii,jj,RemainingChannels(cc)}.ChannelID  = RemainingChannels(cc);
                ops1{ii,jj,RemainingChannels(cc)} = ops;
                ops1{ii,jj,RemainingChannels(cc)}.Ly   = Ly;
                ops1{ii,jj,RemainingChannels(cc)}.Lx   = Lx;
                % use ops.dsprealign to align other images
            end
        end
    end
    
end
%%
fid = cell(numPlanes, ops.NumSplitViews,length(Ca_channel));
for cc=1:ops.nchannels
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
%         fprintf('Copying %s to local (%s)\n',fs{kk}(ff).name,ops.temp_tiff);
%         copyfile(fs{kk}(ff).name,ops.temp_tiff);
        
        ops.temp_tiff = fs{kk}(ff).name;
        
        nFr = nFramesKH(ops.temp_tiff);
        if mod(nFr, nchannels) ~= 0
            fprintf('WARNING: number of frames in tiff (%d) is NOT a multiple of number of channels!\n', ff);
        end
        nChunk = 6000*nchannels;
        nSteps = ceil(nFr/nChunk);
        StartFrames =1+[0:nChunk:nChunk*(nSteps-1)];
        EndFrames   = [nChunk:nChunk:nChunk*nSteps];
        EndFrames(end)=min(nFr,EndFrames(end));
        
        for ss=1:nSteps
            % Because the data is read as a chunk, starting position is not
            % always plane 1. 
            
            for cc=[Ref_channel,RemainingChannels] % start from Ref_channel, after making ops1.dsall
                if ops.doRegistration
                 
                    scan_index = (StartFrames(ss)+cc-1):nchannels:EndFrames(ss);
  
                    % Calculate Reference Image for each channel, based on the movement in reference channel.
                    
                    % Reference Image: mean registered image of its own channel
                    
                    % reference channel: channel within in which movements are used as
                    % a reference for other channels. 
                    
                    if cc==Ref_channel
                        [ops1,data,dsall] = partial_reg(ops,ops1,cc,scan_index,iplane0);
                        subplot(122);plot(dsall(:,1),dsall(:,2),'.-');
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
                            switch ops.RegShape
                                case 'Square'
                                    dreg(yFOVs(:,ll), xFOVs(:,ll), indxr)        = ...
                                        register_movie(data(yFOVs(:,ll), xFOVs(:,ll), indxr), ...
                                        ops1{1,ll,cc}, dsall(indxr,:,ll));
                                case {'GrinLens','TopLeft','TopCenter','TopRight','CenterLeft'}
                                    dreg(:,:, indxr)        = ...
                                        register_movie(data(:,:, indxr), ...
                                        ops1{1,ll,cc}, dsall(indxr,:,ll));
                                otherwise
                                    error('Unknown regshape option %s',ops.RegShape);
                            end
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
                    indframes = ifr0:nplanes:size(dreg,3);
                    for ll = 1:ops.NumSplitViews
                        switch ops.RegShape
                            case 'Square'
                                dwrite = dreg(yFOVs(:,ll),xFOVs(:,ll),indframes);
                            case {'GrinLens','TopLeft','TopCenter','TopRight','CenterLeft'}
                                dwrite = dreg(:,:,indframes);
                            otherwise
                                error('Unknown regshape option %s',ops.RegShape);
                        end
                       
                        datatype = class(dreg);
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
%         ops1{ii}
        if ~any(ops1{ii}.ChannelID == Ca_channel)
            continue; % this channel is not Ca-image. Don't bother to save the registered data. 
        else
            
            fprintf('Loading %s...\n',ops1{ii}.RegFile);
            fid{ii}           = fopen(ops1{ii}.RegFile, 'r');
            
            if ~isempty(ops.RegFileTiffLocation)
                nFramesWrite = 2000;
                ops1{ii} = write_reg_to_tiff2(fid{ii}, ops1{ii}, ii,datatype,nFramesWrite);
                nFramesWrite= 1000; % avoid weird OUT OF MEMORY error.
                write_reg_tave_to_tiff2(fid{ii}, ops1{ii}, ii,ops.RegFile_xtime,datatype,nFramesWrite);
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
    % added to use sorcery
     for ii = 1:numel(ops1)
         ops1{ii}.badframes = false(1, size(ops1{ii}.DS,1));
         if isfield(ops, 'badframes0') && ~isempty(ops.badframes0)
             ops1{ii}.badframes(ops.badframes0) = true;
         end
     end
        
    %%
    % compute xrange, yrange
    for ii = 1:numel(ops1)
        if ops.doRegistration && (ops1{ii}.IsRefChannel)
            minDs = min(ops1{ii}.DS(2:end, [1 2]), [], 1);
            maxDs = max(ops1{ii}.DS(2:end, [1 2]), [], 1);
            disp([minDs(1) maxDs(1) minDs(2) maxDs(2)]);
            if ops.BiDiPhase>0
                maxDs(2) = max(1+ops.BiDiPhase, maxDs(2));
            elseif ops.BiDiPhase<0
                minDs(2) = min(ops.BiDiPhase, minDs(2));
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
    %%
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
            %  a mechanism to work on a subregion of the image 
            dat = data(ops.yFOVs(:,ll),ops.xFOVs(:,ll),indframes);
            if ~isempty(ops.smooth_time_space)
                dat = smooth_movie(dat, ops);
            end
            %<ToDo>: apply log scaling to the data 
            if (ops.IMGLogScaling), dat = log(single(dat)+1); end
            [ds, Corr]  = registration_offsets_KH(dat, ops1{ii,ll,ch_to_align}, 0,ops1{ii,ll,ch_to_align}.RefImg);
            dsall(indframes,:, ll)  = ds;
%             % collect ds; wait, I think it is not necessary. 
%             if ff==1
%                 ds(1,:,:) = 0;
%             end
            ops1{ii,ll,ch_to_align}.DS          = cat(1, ops1{ii,ll,ch_to_align}.DS, ds);
            ops1{ii,ll,ch_to_align}.CorrFrame   = cat(1, ops1{ii,ll,ch_to_align}.CorrFrame, Corr);
        end
   
        
        % added by KH 20161122 to suppress extremely large movement.
        MovementShift = sqrt(sum(ops1{ii,ll,ch_to_align}.DS.^2,2));
        TooBigMovementIndex = find(MovementShift>ops.MaxMovementPixel);
        for jj=1:length(TooBigMovementIndex)
            if TooBigMovementIndex(jj)>1 && TooBigMovementIndex(jj)<length(ops1{ii,ll,ch_to_align}.DS)
                fprintf('!');
                ops1{ii,ll,ch_to_align}.DS(TooBigMovementIndex(jj))= 0.5*( ops1{ii,ll,ch_to_align}.DS(TooBigMovementIndex(jj)-1) +  ops1{ii,ll,ch_to_align}.DS(TooBigMovementIndex(jj)+1));
            end
        end
    end
