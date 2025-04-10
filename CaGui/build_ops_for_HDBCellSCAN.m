function ops = build_ops_for_HDBCellSCAN(db, ops)

ops.build_ver= mfilename;
% ops = db;
ops = addfields(ops, db);
if  ischar(db.expts)
        expts=db.expts;
        db=rmfield(db,'expts')
        db.expts{1}   = expts;
end

% file structure is animalID\sessionID\expts(=SubDirs)\
for k = 1:length(db.expts)
    if iscell(db.expts)
        ops.SubDirs{k}   = db.expts{k};
    elseif isvector(db.expts)
        ops.SubDirs{k}    = num2str(db.expts(k));
    end
end

ops.RootDir = fullfile(ops.RootStorage, ops.mouse_name, ops.date);

Do_Convert_OIR_To_Tiff = 0;

% build file list
ops.fsroot = [];
if isempty(db.RegDirs) % use the file under ops.RootDir\ops.SubDirs\
    for j = 1:length(ops.SubDirs)
        ops.fsroot{j} = dir(fullfile(ops.RootDir, ops.SubDirs{j}, '*.tif'));
        if isempty(ops.fsroot{j})
            if ~exist(fullfile(ops.RootDir, ops.SubDirs{j}),'dir')
                fprintf('No such directory found (%s)',fullfile(ops.RootDir,ops.SubDirs{j}));
            else
                    warning('No .tif files in %s',fullfile(ops.RootDir,ops.SubDirs{j}));
                    fprintf('Found OIR file: Try reading it..,\n');

                oirfile = dir(fullfile(ops.RootDir, ops.SubDirs{j}, '*.oir'));
                if length(oirfile)>1,                    
                    error('There are more than one oir file in %s',fullfile(ops.RootDir, ops.SubDirs{j})); end
                if ~isempty(oirfile) & (Do_Convert_OIR_To_Tiff) %somehow, 32GB OIR file took more than 120GB memory, and PC freezes. Stopping for now.
                    out = ReadImage6D(fullfile(oirfile.folder,oirfile.name), false);
                    image6d=uint16(out{1});
                    image6d=squeeze(permute(image6d,[5,6,4,3,2,1]));
%                     metadata=out{2};
                    [~,TiffName,ext]=fileparts(oirfile.name);
                    TiffName=[TiffName,'.tif'];
                    fprintf('Converting from OIR to Tif file:%s\n',fullfile(oirfile.folder,TiffName));
                    TiffWriter2(image6d,fullfile(oirfile.folder,TiffName),16,'w8'); % bigtiff option
                    clear image6d
                end
            end
        end
        for k = 1:length(ops.fsroot{j})
            ops.fsroot{j}(k).name = fullfile(ops.RootDir, ops.SubDirs{j}, ops.fsroot{j}(k).name);
        end
    end
else
    disp('Use previously registered file');
     for j = 1:length(ops.SubDirs)
        ops.fsroot{j} = dir(fullfile(ops.RootDir, ops.SubDirs{j}, ops.RegDirs{j},'*.tif'));
        if isempty(ops.fsroot{j})
            warning('No .tif files in %s',fullfile(ops.RootDir,ops.SubDirs{j}));
        end
        for k = 1:length(ops.fsroot{j})
            ops.fsroot{j}(k).name = fullfile(ops.RootDir, ops.SubDirs{j}, ops.RegDirs{j}, ops.fsroot{j}(k).name);
        end
    end
end


try    
     % MK code for automatically determining number of planes and channels
    [~, header] = loadFramesBuff(ops.fsroot{1}(1).name, 1, 1, 1);
    
    hh=header{1};
    str = hh(strfind(hh, 'channelsSave = '):end);
    ind = strfind(str, 'scanimage');
    ch = str2num(str(16 : ind(1)-1));
    ops.nchannels = length(ch);
    
    fastZEnable = sscanf(hh(findstr(hh, 'fastZEnable = '):end), 'fastZEnable = %d');
    fastZDiscardFlybackFrames = sscanf(hh(findstr(hh, 'fastZDiscardFlybackFrames = '):end), 'fastZDiscardFlybackFrames = %d');
    if isempty(fastZDiscardFlybackFrames)
        fastZDiscardFlybackFrames = 0;
    end
    stackNumSlices = sscanf(hh(findstr(hh, 'stackNumSlices = '):end), 'stackNumSlices = %d');
    
    ops.nplanes = 1;
    if fastZEnable
        ops.nplanes = stackNumSlices+fastZDiscardFlybackFrames;
    end
    
    str = hh(strfind(hh, 'scanZoomFactor = '):end);
    ind = strfind(str, 'scanimage');
    ops.zoomMicro = str2double(str(18 : ind(1)-1));
catch
end

if ~(isfield(ops, 'planesToProcess') && ~isempty(ops.planesToProcess))
    ops.planesToProcess = 1:ops.nplanes;
else
    % planesToProcess is not working right now
    ops.planesToProcess = 1:ops.nplanes; 
end

CharSubDirs = '';
for i = 1:length(ops.SubDirs)
    CharSubDirs = [CharSubDirs ops.SubDirs{i} '_'];
end
CharSubDirs = CharSubDirs(1:end-1);

ops.ResultsSavePath = sprintf('%s//%s//%s//%s//', ops.ResultsSavePath, ops.mouse_name, ops.date, ...
        CharSubDirs);
    