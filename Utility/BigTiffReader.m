function [F,info] = BigTiffReader(f,index,ind_i,ind_j,ShowProgressBar)
% Tiff reader that can read more than 4GB Tiff file.
% Usage:  [F,info] = BigTiffReader(f)
% 
% --- input ----
% f:  filename 
% 
% --- output ---
% F: tiff data 
% info: information
% 
% Usage2:  [F,info] = BigTiffReader(f,index)
% index: index to read (but you need to know the last index).
% ex): [F,info] = BigTiffReader(f,1:10:100)
% 
% Usage3:  [F,info] = BigTiffReader(f,index,ind_y,ind_x)
% rectangular region indexed by ind_y, ind_x
% ex) rectangular region of F(100:110, 120:150). 
% >> [F,info] = BigTiffReader(f,[],100:110, 120:150)
% 
% 
% Usage4: hide progress bar.
% ShowProgressBar = 0;
% [F,info] = BigTiffReader(f,[],[],[],ShowProgressBar)
% 
% 
% 
% by Kosuke Hamaguchi 20170406

if nargin==0
    f = [];
end


if isempty(f)
    [filename,pathname]=uigetfile({'*.tif'},'Pleae pick a Tiff file');
    f = fullfile(pathname,filename);
else
    [pathname,filename,ext]=fileparts(f);
end

if ~exist(f,'file')
    error('%s not found!',f);
end


info = imfinfo(f);
W=info(1).Width;
H=info(1).Height;


if length(info)<2
    stripOffset = info(1).StripOffsets;
    stripByteCounts = info(1).StripByteCounts;
    Nframes=floor(info(1).FileSize/stripByteCounts);
   
    stripOffset = stripOffset(1) + (0:1:(Nframes-1))'.*stripByteCounts+1;
    stripByteCounts = stripByteCounts*ones(Nframes,1);
  
else
    Nframes=length(info);
    stripOffset = cat(1,info.StripOffsets);
% stripOffset = info(1).StripOffsets;
    stripByteCounts = cat(1,info.StripByteCounts);
%     stripByteCounts = stripByteCounts(:,1)';
end

  RowsPerStrip = info(1).RowsPerStrip;

if nargin<=1
    index = 1:Nframes;
elseif isempty(index)
    index = 1:Nframes;
end

if max(index)>Nframes
    warning('index(end)=%d > Nframes(%d). index(end) is limited to data size %d',...
        index(end), Nframes,Nframes);
    index= index(index<=Nframes);
end

if nargin<=2 || isempty(ind_i),     ind_i = 1:H;    end    
if nargin<=3  || isempty(ind_j),    ind_j = 1:W;    end
if nargin<=4  || isempty(ShowProgressBar),    ShowProgressBar = 1;    end

if info(1).BitDepth==32
    typedef = 'uint32=>uint32';
elseif info(1).BitDepth==16
    typedef = 'uint16=>uint16';
elseif info(1).BitDepth==8
    typedef =  'uint8=>uint8';
else
    error('Unknown bit depth %d',info(1).BitDepth)
end


%%

switch info(1).ByteOrder
    case  'big-endian'
        fID = fopen (f, 'r','b');
        read_pos_offset = 0;
    case 'little-endian'
        fID = fopen (f, 'r','l');
        read_pos_offset = 0;
    otherwise
        error('Unknown ByteOrder %s',info(1).ByteOrder);
end

switch info(1).BitDepth
    case 32
        F = zeros(length(ind_i),length(ind_j),length(index),'uint32');
        tmp =  uint32(zeros(W,H)); % data is stored in transposed view
    case 16
        F  = zeros(length(ind_i),length(ind_j),length(index),'uint16');
        tmp =  uint16(zeros(W,H)); % data is stored in transposed view
    case 8
        F = zeros(length(ind_i),length(ind_j),length(index),'uint8');
        tmp =  uint8(zeros(W,H)); % data is stored in trasposed view
end

if ShowProgressBar
waitCnt = 0;
waitH=waitbar(waitCnt,'Loading Tiffs ( / )...');
end

for ii = 1:length(index)
    
    if ShowProgressBar &&  ii/length(index)>=waitCnt
        waitbar(waitCnt,waitH,{strrep(filename,'_','\_'),sprintf('Loading Tiffs (%d/%d)',...
            ii,length(index))});
        waitCnt=waitCnt+0.1;
        %nMsgChars = overfprintf(nMsgChars, '%i/%i', t, nFrames);
    end
    
%         read_pos = stripOffset(index(ii),:)-info(1).Offset/8;
         read_pos = stripOffset(index(ii),:)+read_pos_offset;
%     read_pos = stripOffset(index(ii),:);
    read_bytes = stripByteCounts(index(ii),:);
    for jj=1:length(read_pos)
%         H = info(index(ii)).RowsPerStrip;
        W = read_bytes(jj)/(info(1).BitDepth/8*RowsPerStrip);
        rows_index = [1:RowsPerStrip]+(jj-1)*RowsPerStrip;
        fseek(fID, read_pos(jj), 'bof');
       
%         tmp(rows_index,:) = fread(fID, [H,W],typedef);
         tmp(:,rows_index) = fread(fID, [W,RowsPerStrip],typedef);
        
         
    end

        F(:,:,ii)=tmp(ind_j,ind_i)';

end

 fclose(fID);
if ShowProgressBar, close(waitH); end