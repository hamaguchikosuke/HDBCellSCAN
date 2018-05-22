function IMG = get_sparse_subset_images(SubDirs,fs,nchannels,nplanes,ChannelSelected,ops)

BiDiPhase            = getOr(ops, {'BiDiPhase'}, 0);

nbytes = fs{1}(1).bytes;
nFr = nFramesKH(fs{1}(1).name); % assuming that all the tiffs have the same nunber of images

ntifs = sum(cellfun(@(x) numel(x), fs));
nfmax = min(ceil(ops.NimgFirstRegistration/ntifs),1999);

indx = 0;
if nFr<(nchannels*nplanes*nfmax)
    warning('NimgFirstRegistration(%d) is larger than number of images per plane per channel');
    ops.NimgFirstRegistration = floor(ntifs*nFr/nplanes/nchannels);
    fprintf('ops.NimgFirstRegistration is set to %d\n',ops.NimgFirstRegistration);
end

        
IMG = zeros(ops.Ly, ops.Lx, nplanes, ops.NimgFirstRegistration, 'int16');

for k = 1:length(SubDirs)
    iplane0 = 1;
    for j = 1:length(fs{k})
        if abs(nbytes - fs{k}(j).bytes)>1e3
            nbytes = fs{k}(j).bytes;
            nFr = nFramesKH(fs{k}(j).name); % nFrame is the function to obtain number of frames in Tiff.
        end
       
        iplane0 = mod(iplane0-1, nplanes) + 1;
%         offset = 0;
        if j==1
            offset=0;
%             offset = nchannels*(nplanes-1);
        end
        ichanset = [offset + nchannels*(iplane0-1) + [ChannelSelected;...
            nchannels*nplanes*nfmax]; nchannels];
        
        
        % ichanset: [START,END,SKIP]
        fprintf('Loading frame [%d:%d:%d]\n',ichanset(1),ichanset(3),ichanset(2));
        %             data = loadFramesBuff2(fs{k}(j).name, ichanset(1),ichanset(2), ichanset(3));
        data = BigTiffReader(fs{k}(j).name, [ichanset(1):ichanset(3):ichanset(2)]);
        data = reshape(data, size(data,1), size(data,2), nplanes, []);
        
        if BiDiPhase % to correct bi-directional phase shifting
            data = ResonanceImagingPhaseShifter(data,BiDiPhase);
        end
        IMG(:,:,:,indx+(1:size(data,4))) = data;
        indx = indx + size(data,4);
        
        iplane0 = iplane0 - nFr/nchannels;
    end
end
IMG =  IMG(:,:,:,1:indx);