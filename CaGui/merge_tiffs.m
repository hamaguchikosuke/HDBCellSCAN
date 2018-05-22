% merge all tiffs 
TiffPath=uigetdir();
S=dir(fullfile(TiffPath,'*.tif'));
[~,sort_ind]=sort([S.datenum]);
S=S(sort_ind);

 nFramesWrite = 2000;
 WriteCnt = nFramesWrite:nFramesWrite:length(S);
 tmp = imread(fullfile(TiffPath,S(1).name));
 switch class(tmp)
     case 'uint8'
         data= uint8(zeros(size(tmp,1),size(tmp,2),nFramesWrite));
          bitspersamp = 8;
     case 'uint16'
         data= uint16(zeros(size(tmp,1),size(tmp,2),nFramesWrite));
          bitspersamp = 16;
 end
 cnt=1;
 fname = [strrep(S(1).name,'.tif',''),'_merged.tif'];
 write_name = fullfile(TiffPath,fname);
 for ii=1:length(S)
     tmp = imread(fullfile(TiffPath,S(ii).name));
     data(:,:,cnt)=tmp;
     cnt=cnt+1;
     if any(WriteCnt==ii)
         TiffWriter(data,write_name,bitspersamp);
         cnt=1;
         fprintf('%d/%d...',ii,length(S));
     end
 end
 
 % Write the final chunk of tiff
 if cnt~=1
     data=data(:,:,1:cnt-1);
      TiffWriter(data,write_name,bitspersamp);
 end
 fprintf('Done\n');