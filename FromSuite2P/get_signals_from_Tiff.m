function Fcell = get_signals_from_Tiff(ops, iplane)
loadname = sprintf('%s/F_%s_%s_plane%d_Nk%d.mat',...
       ops.ResultsSavePath,ops.mouse_name, ops.date, iplane, ops.Nk);
try
   load(loadname); % load ops, res, stat,stat0,res0,clustrules.
catch
   error('Could not find cell detection file %s\n',loadname); 
end

TiffFiles=[];
EstimNumSeries = [];
for ii=1:length(ops.SubDirs)
    RegFastFilePath=fullfile(ops.RegFileTiffLocation,ops.mouse_name,ops.date,ops.SubDirs{ii},sprintf('Plane%d',ops.iplane));
    RegFastFile=sprintf('%s_%s_%s_2P_plane%d_*.tif',ops.date, ops.SubDirs{ii},ops.mouse_name,ops.iplane);
      f = dir(fullfile(RegFastFilePath,RegFastFile));
      
    if ~isempty(f)
       ftmp = fullfile(RegFastFilePath,{f.name});
    else
        error('reg files not found: %s',fullfile(RegFastFilePath,RegFastFile));
    end
    TiffFiles = cat(1,TiffFiles,ftmp);
    tt = Tiff(TiffFiles{1});
    w=tt.getTag('ImageWidth');
    h=tt.getTag('ImageLength');
    EstimNumSeries = cat(1,EstimNumSeries,floor([f.bytes]/tt.getTag('BitsPerSample')*8/(w*h)));
end
CumNSeries = cumsum(EstimNumSeries);

if ops.Nframes(end)~=CumNSeries(end)
    error('Total number of frame does not match');
end
Nk = numel(stat);
% Nkpar = ops.Nk;


%% get signals  
[Ly Lx] = size(ops.mimg);

% nimgbatch = 2000;

ind_start = 1+[0 CumNSeries(1:end-1)];
ind_end =    CumNSeries;


tic
F = zeros(Nk, sum(ops.Nframes), 'single');
cnt = 0;
for ii=1:length(TiffFiles)
   
    data = single(loadFramesBuff2(TiffFiles{ii}));
    data = data(ops.yrange, ops.xrange, :);
    NT= size(data,3);
    data = reshape(data, [], NT);
    index = cnt+(1:NT);
    for k = 1:Nk
        ipix = stat(k).ipix;
        if ~isempty(ipix)
            %            F(k,ix + (1:NT)) = stat(k).lambda' * data(ipix,:);
            F(k,index) = mean(data(ipix,:), 1);
        end
    end
    
    fprintf('Frame %d/%d done in time %2.2f \n', ii,length(TiffFiles), toc)
    cnt = index(end);
end

csumNframes = [0 cumsum(ops.Nframes)];
Fcell = cell(1, length(ops.Nframes));
for i = 1:length(ops.Nframes)
    Fcell{i} = F(:, csumNframes(i) + (1:ops.Nframes(i)));
end

save(sprintf('%s/F_%s_%s_plane%d_Nk%d.mat', ops.ResultsSavePath, ...
    ops.mouse_name, ops.date, iplane, ops.Nk),  'ops', 'res', 'stat', 'stat0', 'res0', 'Fcell', 'clustrules')

