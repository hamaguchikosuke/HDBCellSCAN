function out = TaveTiffs(f,nAve)
% Write N-frame averaged Tiff. 
% Usage:  TaveTiffs(f,nAve);
% f: full filename or empty
% nAve: N frames to be averagted.
% 
% ex)  TaveTiffs([],5) % 5-frame average 
% 
%  by KH 20170406

if isempty(f)
[f,TgtPath]=uigetfile(fullfile(pwd,'*.tif'),'Pleae select target Tiffs','MultiSelect','on');
else
    [TgtPath,f,ext]=fileparts(f);
end

SavePath = fullfile(TgtPath,sprintf('x%dmovie',nAve));
if ~exist(SavePath,'dir')
mkdir(SavePath);
end

if isempty(f)
    error('No .tif file defined!');
end

if ischar(f)
    f = {f};
end

SaveCandidate = f{1};
[SaveName,SavePath]=uiputfile(fullfile(SavePath,SaveCandidate),'Save time-averaged Tiff as (it will add x5_001.tif)');


% <To Do>
% Make sure the file order is correct.

% RegFileSearchString=fullfile(TgtPath,'*.tif');
% f=dir(RegFileSearchString);


Nframes_Each = zeros(1,length(f));

for ff=1:length(f)
    
    FullTiffName=fullfile(TgtPath,  f{ff});
    
    info = imfinfo(FullTiffName);
    W=info(1).Width;
    H=info(1).Height;
    
    try
        pat = '\w*frames=(?<nFrames>\d+)\n\w';
        Nframes=regexp(info.ImageDescription,pat,'tokens');
        Nframes_Each(ff) =   str2double(Nframes{1});
    catch
        Nframes_Each(ff) = length(info);
    end
end

NTotalFrames = sum(Nframes_Each);
fprintf('In total %d frames\n',NTotalFrames)

Nbatch = 6000;
NTiffWrite = 2000;

if mod(NTotalFrames,nAve)~=0
    error('Number of total frames %d is not divisive of nAve %d',...
        NTotalFrames,nAve);
end
%%

fprintf('===Writing x%d movie ====\n',nAve);    
    TiffData = uint16(zeros(H,W,NTiffWrite));
    TiffIndex= 0;
    SaveCnt = 0;
    
for ff=1:length(f)
    
    FullTiffName=fullfile(TgtPath,  f{ff});
    
    info = imfinfo(FullTiffName);
    BitsPerSampe = info.BitsPerSample;
%     W=info1(1).Width;
%     H=info1(1).Height;
     
    Nframes =  Nframes_Each(ff) ;
    
    NStride = ceil(Nframes/Nbatch);
    BatchTiff_Start = 1+[0:Nbatch:(NStride-1)*Nbatch];

  
    
    for ii = 1:NStride
        index = BatchTiff_Start(ii)+[0:Nbatch-1];
        index = index(index<=Nframes);
        fprintf ('\nload image ... %d...%d,', index(1), index(end));
        F = BigTiffReader(FullTiffName,index);
        
        [H,W,T]=size(F);
        Fz = reshape(F,[H,W,nAve,T/nAve]); 
        Fz = squeeze(mean(Fz,3));
        
        
        
        TiffIndex = TiffIndex(end) + [1:size(Fz,3)];
        
        if TiffIndex(end)>NTiffWrite % in case it does not fit 
            remain_ind= 1:size(Fz,3);
            remain_ind = remain_ind(TiffIndex>NTiffWrite);
            save_ind =   1:size(Fz,3);
            save_ind =   save_ind(TiffIndex<=NTiffWrite);
            Fz_remain = Fz(:,:,remain_ind);
            Fz        = Fz(:,:,save_ind);
            TiffIndex = TiffIndex(TiffIndex<=NTiffWrite);
        else
            Fz_remain = [];
        end
        
        TiffData(:,:,TiffIndex)=Fz;
        TiffIndex = TiffIndex(end);
  
        
        if TiffIndex(end)==NTiffWrite
            fprintf(' Saving Tiff...');
            SaveCnt = SaveCnt+1;
     
            [~,TiffName,Ext]=fileparts(SaveName);
            ExportTiffName = sprintf('%s_x%d_%03d.tif',TiffName,nAve,SaveCnt);
            fname = fullfile(SavePath, ExportTiffName);
            
            TiffWriter(uint16(TiffData),fname,BitsPerSampe);
            TiffIndex = 0;
        end
        
        if ~isempty(Fz_remain)
            TiffIndex = TiffIndex(end) + [1:size(Fz_remain,3)];
            TiffData(:,:,TiffIndex)=Fz_remain;
            TiffIndex = TiffIndex(end);
        end
        
    end
    
   
    
end

if TiffIndex(end)~=NTiffWrite & TiffIndex(end)>0 % if end with exact number, TiffIndex(end) is zero.
    fprintf('Saving Remaining Tiff...\n');
    SaveCnt = SaveCnt+1;
    [~,TiffName,Ext]=fileparts(SaveName);
    ExportTiffName = sprintf('%s_x%d_%03d.tif',TiffName,nAve,SaveCnt);
    fname = fullfile(SavePath, ExportTiffName);
    TiffData = TiffData(:,:,1:TiffIndex(end));
    TiffWriter(uint16(TiffData),fname,BitsPerSampe);
end

fprintf('====  Done =====\n')

out = 1;
end
