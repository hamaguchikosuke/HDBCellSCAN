function ops = write_reg_tave_to_tiff2(fid, ops, iplane,xtime,datatype,varargin)
% write time-averaged, registered (motion-corrected) tiff files 
% xtime: an integer how many frames to combine.

if nargin>=6
    nFramesWrite = varargin{1};
else
    nFramesWrite = 2000;
end

Ly = ops.Ly;
Lx = ops.Ly;
bitspersamp = 16;

frewind(fid);
for k = 1:length(ops.SubDirs)
    ix = 0;    
    nframesleft = ops.Nframes(k);
%     if mod(ops.Nframes(k),xtime)==0
%         % the number of tiffs are multiple integeres of xtime, Good.
%     else
%         error('the number of tiffs %d are NOT multiple integeres of xtime = %d',ops.Nframes(k),xtime);
%     end
    datend = [];
    while nframesleft>0
         ix = ix + 1;
         NFrames2Read = ceil(nFramesWrite/xtime)*xtime*5; % roughly 10000 frames per read.
        
         if nframesleft >= NFrames2Read % if there are enough tiffs left
            nfrtoread = min(nframesleft, NFrames2Read);
         else % if remaining tiffs are less than read, limit the number of read to the multiple integer of xtime.
            nfrtoread = xtime*floor(nframesleft/xtime);
         end
         if nfrtoread==0
             break;
         end
        data = fread(fid,  Ly*Lx*nfrtoread, datatype);                
        nframesleft = nframesleft - nfrtoread;
       
        data = reshape(data, Ly*Lx, xtime, []);        
        data = mean(data,2);
        data = reshape(data,Ly,Lx,[]);
        
        [~,subfolder,~]=fileparts(ops.RegFile);
        subfolder = strrep(subfolder,'tempreg_','');
        foldr = fullfile(ops.ResultsSavePath, sprintf('%s\\x%dmovie',subfolder,xtime));
        

        if ~exist(foldr, 'dir')
            mkdir(foldr)
        end
        partname = sprintf('%s_%s_%s_2P_%s_x%d_%03d.tif', ops.date, ops.SubDirs{k}, ...
            ops.mouse_name, subfolder, xtime,ix);
        fname = fullfile(foldr, partname);
        
        writeopt = 'w'; % append
        fprintf('Writing time averaged Tiff %d\n',ix);
        TiffWriter2(uint16(data),fname,bitspersamp,writeopt);
    end
end
