function folder = get_tiff_location_from_ops(ops,reg2P_ver)
% return registered movie folder from ops.
% folder = get_tiff_location_from_ops(ops,reg2P_ver);
% ops: options 
% reg2P_ver: name of the reg2P program used to generate registered movie.

switch reg2P_ver
    case {'reg2P_kh004','reg2P_kh005'}
        [~,subfolder,~]=fileparts(ops.RegFile);
        subfolder = strrep(subfolder,'tempreg_','');
        folder = fullfile(ops.ResultsSavePath,subfolder);
    case 'reg2P_kh003'
        error('not defined yet!');
    otherwise
        error('Unknown reg2P version %s~',reg2P_ver)
end