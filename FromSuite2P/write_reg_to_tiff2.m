function ops = write_reg_to_tiff2(fid, ops, iplane,datatype, varargin)

if nargin>=5
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
    
    datend = [];
    while nframesleft>0
         ix = ix + 1;
        nfrtoread = min(nframesleft, nFramesWrite);
        data = fread(fid,  Ly*Lx*nfrtoread, datatype);                
        nframesleft = nframesleft - nfrtoread;
        data = reshape(data, Ly, Lx, []);        

        [~,subfolder,~]=fileparts(ops.RegFile);
        subfolder = strrep(subfolder,'tempreg_','');
        foldr = fullfile(ops.RegFileTiffLocation, ops.mouse_name, ops.date, ...
            ops.SubDirs{k}, subfolder);
        if ~exist(foldr, 'dir')
            mkdir(foldr)
        end
        partname = sprintf('%s_%s_%s_2P_%s_%03d.tif', ops.date, ops.SubDirs{k}, ...
            ops.mouse_name, subfolder, ix);
        fname = fullfile(foldr, partname);
        
        TiffWriter(uint16(data),fname,bitspersamp);
    end
end
