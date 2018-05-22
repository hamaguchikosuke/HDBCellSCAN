function ops1 = reg2P_kh002(ops)
%%
%  removed ichannel and rchannel, as we use RCaMP and ERK. 
numPlanes = length(ops.planesToProcess);

nplanes             = getOr(ops, {'nplanes'}, 1);
nchannels           = getOr(ops, {'nchannels'}, 1);

Ca_channel          = getOr(ops, {'Ca_channel'}, 1);
Align_channel         = getOr(ops, {'Align_channel'}, 2);    

% red_align           = getOr(ops, {'AlignToRedChannel'}, 0);
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

nbytes = fs{1}(1).bytes;
nFr = nFrames(fs{1}(1).name);
Info0 = imfinfo(fs{1}(1).name);
Ly = Info0(1).Height;
Lx = Info0(1).Width;
ops.Lx = Lx;
ops.Ly = Ly;

[xFOVs, yFOVs] = get_xyFOVs(ops); 

indx = 0;
IMG = zeros(Ly, Lx, nplanes, ops.NimgFirstRegistration, 'int16');

if ops.doRegistration
    for k = 1:length(ops.SubDirs)
        iplane0 = 1;
        for j = 1:length(fs{k})
            if abs(nbytes - fs{k}(j).bytes)>1e3
                nbytes = fs{k}(j).bytes;
                nFr = nFrames(fs{k}(j).name); % nFrame is the function to obtain number of frames in Tiff.
            end
            if nFr<(nchannels*nplanes*nfmax + nchannels*nplanes)
                continue;
            end
            
            iplane0 = mod(iplane0-1, nplanes) + 1;
            offset = 0;
            if j==1
                offset = nchannels*nplanes;
            end
            ichanset = [offset + nchannels*(iplane0-1) + [Align_channel;...
                nchannels*nplanes*nfmax]; nchannels];
            
            iplane0 = iplane0 - nFr/nchannels;
            % ichanset: [START,END,SKIP]
            fprintf('Loading frame [%d:%d:%d]\n',ichanset(1),ichanset(3),ichanset(2));
            data = BigTiffReader(fs{k}(j).name, [ichanset(1):ichanset(3):ichanset(2)]);
            data = reshape(data, Ly, Lx, nplanes, []);
            
            if BiDiPhase % to correct bi-directional phase shifting
                data = ResonanceImagingPhaseShifter(data,BiDiPhase);
            end
            IMG(:,:,:,indx+(1:size(data,4))) = data;
            indx = indx + size(data,4);
            
        end
    end
    IMG =  IMG(:,:,:,1:indx);
    
    ops1 = cell(numPlanes, size(xFOVs,2));
    for i = 1:numPlanes
        for j = 1:size(xFOVs,2)
            ops1{i,j} = align_iterative(single(squeeze(IMG(yFOVs(:,j),xFOVs(:,j),...
                ops.planesToProcess(i),:))), ops);
        end
    end
    
    if ops.showTargetRegistration
        h=figure(11);
        set(h,'position', [900 50 900 300*numPlanes])
        ax = ceil(sqrt(numel(ops1)/2));
        i0 = 0;
        for i = 1:numPlanes
            for j = 1:size(xFOVs,2)
                i0 = i0+1;
                subplot(ax,2*ax,i0)
                imagesc(ops1{i,j}.mimg)
                colormap('gray')
                title(sprintf('Registration for plane %d, mouse %s, date %s', ...
                    i, ops.mouse_name, ops.date))
            end
        end
        
        drawnow
    end
    
    clear IMG
else
    for i = 1:numPlanes
        ops1{i} = ops;
        ops1{i}.mimg = zeros(Ly, Lx);
        ops1{i}.Ly   = Ly;
        ops1{i}.Lx   = Lx;
     end
end
%%
fid = cell(numPlanes, size(xFOVs,2));
for i = 1:numPlanes
    for j = 1:size(xFOVs,2)
        ops1{i,j}.RegFile = fullfile(ops.RegFileRoot, sprintf('tempreg_plane%d.bin', i + (j-1)*numPlanes));
        regdir = fileparts(ops1{i,j}.RegFile);
        if ~exist(regdir, 'dir')
            mkdir(regdir);
        end
        
        % open bin file for writing
        fid{i,j}              = fopen(ops1{i,j}.RegFile, 'w');
        ops1{i,j}.DS          = [];
        ops1{i,j}.CorrFrame   = [];
        ops1{i,j}.mimg1       = zeros(ops1{i,j}.Ly, ops1{i,j}.Lx);
        
    end
end


%% Then by using the reference frame, align all the frames.
tic
for k = 1:length(fs)
    for i = 1:numel(ops1)
         ops1{i}.Nframes(k)  = 0;
    end
    
    iplane0 = 1:1:ops.nplanes;
    for j = 1:length(fs{k})
        iplane0 = mod(iplane0-1, numPlanes) + 1;
        fprintf('Copying %s to local (%s)\n',fs{k}(j).name,ops.temp_tiff);
        copyfile(fs{k}(j).name,ops.temp_tiff)
        
        nFr = nFrames(ops.temp_tiff);
        if mod(nFr, nchannels) ~= 0
            fprintf('  WARNING: number of frames in tiff (%d) is NOT a multiple of number of channels!\n', j);
        end
        
        nChunk = 6000;
        if rem(nChunk,nchannels)~=0, error('nChunk must be a multiple integers of nchannel(%d)',nchannel); end
        nSteps = ceil(nFr/nChunk);
        StartFrames =1+[0:nChunk:nChunk*(nSteps-1)];
        EndFrames   = [nChunk:nChunk:nChunk*nSteps];
        EndFrames(end)=min(nFr,EndFrames(end));
        
        for ss=1:nSteps
            
        scan_index = (StartFrames(ss)+Align_channel-1):nchannels:EndFrames(ss);
       
        data = BigTiffReader(ops.temp_tiff, scan_index);
        
        
        if BiDiPhase
          data = ResonanceImagingPhaseShifter(data,BiDiPhase);
        end
        
        %---------- do registration -------------%
        if ops.doRegistration
            % get the registration offsets
            dsall = zeros(size(data,3), 2, size(xFOVs,2));
            for i = 1:numPlanes
                ifr0 = iplane0(ops.planesToProcess(i));
                indframes = ifr0:nplanes:size(data,3);

                for l = 1:size(xFOVs,2)
                    dat = data(yFOVs(:,l),xFOVs(:,l),indframes);
                    if ~isempty(ops.smooth_time_space)
                        dat = smooth_movie(dat, ops); 
                    end
                    [ds, Corr]  = registration_offsets(dat, ops1{i,l}, 0);
                    dsall(indframes,:, l)  = ds;
                    % collect ds
                    if j==1
                        ds(1,:,:) = 0;
                    end
                    ops1{i,l}.DS          = cat(1, ops1{i,l}.DS, ds);
                    ops1{i,l}.CorrFrame   = cat(1, ops1{i,l}.CorrFrame, Corr);
                end
            end
            
            % added by KH 20161122 to suppress extremely large movement.
            MovementShift = sqrt(sum(ops1{i,l}.DS.^2,2));
            TooBigMovementIndex = find(MovementShift>ops.MaxMovementPixel);
            for ii=1:TooBigMovementIndex
                if TooBigMovementIndex(ii)>1 || TooBigMovementIndex(ii)<length(ops1{i,l}.DS)
                    fprintf('!');
                ops1{i,l}.DS(TooBigMovementIndex(ii))= 0.5*( ops1{i,l}.DS(TooBigMovementIndex(ii-1)) +  ops1{i,l}.DS(TooBigMovementIndex(ii+1)))
                end
            end
                
            
            % if aligning by the red channel, data needs to be reloaded as the
            % green channel
            if Ca_channel ~= Align_channel
                scan_index = (StartFrames(ss)+Ca_channel-1):nchannels:EndFrames(ss);
                data = BigTiffReader(ops.temp_tiff, scan_index);                 
            end
            
            ix0 = 0;
            Nbatch = 1000;
            dreg = zeros(size(data), class(data)); 
%             fprintf('Preparing dreg zeros(%d,%d,%d)\n',size(data));
            
            while ix0<size(data,3)
                indxr = ix0 + (1:Nbatch);
                indxr(indxr>size(data,3)) = [];
%                 fprintf('indxr=[%d...%d]\n',indxr(1),indxr(end))
                for l = 1:size(xFOVs,2)
                    dreg(yFOVs(:,l), xFOVs(:,l), indxr)        = ...
                        register_movie(data(yFOVs(:,l), xFOVs(:,l), indxr), ops1{1,l}, dsall(indxr,:,l));
                end
                ix0 = ix0 + Nbatch;
            end
            
        else   %----------  NO registration -------------%
            dreg = data;
        end
        
        % write dreg to bin file+
        for i = 1:numPlanes
            ifr0 = iplane0(ops.planesToProcess(i));
            indframes = ifr0:nplanes:size(data,3);
            for l = 1:size(xFOVs,2)
                dwrite = dreg(yFOVs(:,l),xFOVs(:,l),indframes);
                datatype = class(data);
                fwrite(fid{i,l}, dwrite, datatype);
                
                ops1{i,l}.Nframes(k) = ops1{i,l}.Nframes(k) + size(dwrite,3);
                ops1{i,l}.mimg1 = ops1{i,l}.mimg1 + sum(dwrite,3);
            end
        end
     
            fprintf('Chunk %d (%d/%d) \n', ss,EndFrames(ss),nFr)

        end
        if rem(j,5)==1
            fprintf('Set %d, tiff %d done in time %2.2f \n', k, j, toc)            
        end
        
        iplane0 = iplane0 - nFr/nchannels;
    end
    
end
for i = 1:numel(ops1)
    ops1{i}.mimg1 = ops1{i}.mimg1/sum(ops1{i}.Nframes);
end
%%
for i = 1:numPlanes    
    fclose(fid{i});
    fprintf('Loading %s...\n',ops1{i}.RegFile);
    fid{i}           = fopen(ops1{i}.RegFile, 'r');
    
    if ~isempty(ops.RegFileTiffLocation)
        ops1{i} = write_reg_to_tiff(fid{i}, ops1{i}, i,datatype);
        write_reg_tave_to_tiff(fid{i}, ops1{i}, i,ops.RegFile_xtime,datatype);
    end    
    
    if ~isempty(ops.nimgbegend) && ops.nimgbegend>0
        ops1{i} = getBlockBegEnd(fid{i}, ops1{i}); % get mean of first and last frames in block (to check for drift)
    end
    
    if ~isempty(ops.RegFileBinLocation)
        folder = fullfile(ops1{i}.RegFileBinLocation, ops1{i}.mouse_name, ...
            ops1{i}.date);
        if ~exist(folder, 'dir')
            mkdir(folder)
        end
        fidCopy = fopen(fullfile(folder, sprintf('plane%d.bin', i)), 'w');
        sz = ops1{i}.Lx * ops1{i}.Ly;
        parts = ceil(sum(ops1{i}.Nframes) / 2000);
        for p = 1:parts
            toRead = 2000;
            if p == parts
                toRead = sum(ops1{i}.Nframes) - 2000 * (parts-1);
            end
            data = fread(fid{i},  sz*toRead, '*int16');
            fwrite(fidCopy, data, class(data));
        end
        fclose(fidCopy);
    end
    fclose(fid{i});
end
%%
% compute xrange, yrange
for i = 1:numel(ops1)
    if ops.doRegistration
        minDs = min(ops1{i}.DS(2:end, [1 2]), [], 1);
        maxDs = max(ops1{i}.DS(2:end, [1 2]), [], 1);
        disp([minDs(1) maxDs(1) minDs(2) maxDs(2)])
        if BiDiPhase>0
            maxDs(2) = max(1+BiDiPhase, maxDs(2));
        elseif BiDiPhase<0
            minDs(2) = min(BiDiPhase, minDs(2));
        end
        
        ops1{i}.yrange = ceil(maxDs(1)):floor(ops1{i}.Ly+minDs(1));
        ops1{i}.xrange = ceil(maxDs(2)):floor(ops1{i}.Lx+minDs(2));
    else
        ops1{i}.yrange = 1:Ly;
        ops1{i}.xrange = 1:Lx;
    end
    
    savepath = sprintf('%s/', ops.ResultsSavePath);
    
    if ~exist(savepath, 'dir')
        mkdir(savepath)
    end
    ops = ops1{i};
    save(sprintf('%s/regops_%s_%s_plane%d.mat', ops.ResultsSavePath, ...
        ops.mouse_name, ops.date, i),  'ops')
end


%save(sprintf('%s/F_%s_%s_plane%d.mat', ops.ResultsSavePath, ...
 %   ops.mouse_name, ops.date, iplane), 'ops')

%%
