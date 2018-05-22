function out = ApplyRegOps2Tiff()
% Apply x-y movement shift to Tiff stacks. 
% To use this function, load movement data from regops*.mat in which ops.DS
% represents the movements in x-y direction.
% Usage: ApplyRegOps2Tiff()
% First, it will ask regops*.mat file.
% Second, it will ask target Tiff files. They needs to be numbered
% correctly so that movement shift is applied to the correct frames.
% 
% by KH 20170406

[regopsFile,pathname]=uigetfile({'regops*.mat'},'Pleae pick a RegOps file');

FullRegOpsFile = fullfile(pathname,regopsFile);
% FullRegOpsFile = 'G:\Kosuke\Data\F\DualLickMice\b6j594\20170405\1\regops_B6j594_20170405_plane1.mat';
load(FullRegOpsFile);

% To use GPU, set 1
regops.useGPU = 1;

% Tiff is processed for every Nbatch frames as a batch, and 
regops.Nbatch = 1000;

% registered multi-page Tiff file is generated for every NTiffWrite frames
regops.NTiffWrite = 3*regops.Nbatch;

% show movie for every show_movie frames.
regops.show_movie = 10;
%%
% [TgtPath]=uigetdir(pathname,'Pleae select the directory that contains target Tiffs');
[f,TgtPath]=uigetfile(fullfile(pathname,'*.tif'),'Pleae select target Tiffs','MultiSelect','on');

% <To Do>
% Make sure the file order is correct.

% RegFileSearchString=fullfile(TgtPath,'*.tif');
% f=dir(RegFileSearchString);

if isempty(f)
    error('No .tif file defined!');
end

if ischar(f)
    f = {f};
end

Nframes_Each = zeros(1,length(f));

for ff=1:length(f)
    
    FullTiffName=fullfile(TgtPath,  f{ff});
    
    info = imfinfo(FullTiffName);
    W=info(1).Width;
    H=info(1).Height;
    
    pat = '\w*frames=(?<nFrames>\d+)\n\w';
    Nframes=regexp(info.ImageDescription,pat,'tokens');
    Nframes_Each(ff) =   str2double(Nframes{1});
    
end
% sanity check
if sum(Nframes_Each)~=length(ops.DS)
    error('The number of frame is different (Tiff: %d), (movement: %d)\n',...
        sum(Nframes_Each),length(ops.DS));
else
    fprintf('==== Image registration (Total: %d frames) ====\n',sum(Nframes_Each));
end

%%

for ff=1:length(f)
    
    FullTiffName=fullfile(TgtPath,  f{ff});
    
    info = imfinfo(FullTiffName);
    BitsPerSampe = info.BitsPerSample;
%     W=info1(1).Width;
%     H=info1(1).Height;
     
    pat = '\w*frames=(?<nFrames>\d+)\n\w';
    Nframes=regexp(info.ImageDescription,pat,'tokens');
    Nframes =   str2double(Nframes{1});
    
    NStride = ceil(Nframes/regops.Nbatch);
    BatchTiff_Start = 1+[0:regops.Nbatch:(NStride-1)*regops.Nbatch];
    
    TiffData = uint16(zeros(H,W,regops.NTiffWrite));
    TiffCnt = 0;
    NSave = regops.NTiffWrite/regops.Nbatch;
    SaveCnt = 0;
    ii_cnt = 0;
    
    for ii = 1:NStride
        index = BatchTiff_Start(ii)+[0:regops.Nbatch-1];
        index = index(index<=Nframes);
        fprintf ('\nregistering image ... %d...%d,', index(1), index(end));
        F = BigTiffReader(FullTiffName,index);
        ds = ops.DS(index,:);
        dreg = register_movie(F, ops, ds);
        
        TiffIndex = TiffCnt+[1:regops.Nbatch];
        TiffData(:,:,TiffIndex)=dreg;
        TiffCnt = max(TiffIndex);
        ii_cnt = ceil(ii_cnt+1);
        
        if ii_cnt >= NSave 
            fprintf(' Saving Tiff...');
            SaveCnt = SaveCnt+1;
            ii_cnt = 0;
            [~,TiffName,Ext]=fileparts(FullTiffName);
            ExportTiffName = sprintf('%s_%03d.tif',TiffName,SaveCnt);
            fname = fullfile(TgtPath, ExportTiffName);
            
            TiffWriter(uint16(TiffData),fname,BitsPerSampe);
            TiffCnt = 0;
        end
     
    end
    
    if ii_cnt ~=0
        fprintf('Saving Remaining Tiff...\n');
        SaveCnt = SaveCnt+1;
        ii_cnt = 0;
        [~,TiffName,Ext]=fileparts(FullTiffName);
        ExportTiffName = sprintf('%s_%03d.tif',TiffName,SaveCnt);
        fname = fullfile(TgtPath, ExportTiffName);
        TiffData = TiffData(:,:,1:TiffCnt);
        TiffWriter(uint16(TiffData),fname,BitsPerSampe);
        TiffCnt = 0;
    end
    
end

fprintf('====  Done =====\n')

end

