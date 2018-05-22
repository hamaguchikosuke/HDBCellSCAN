function [stat, res, ops] = apply_ROIrules_kh(ops, stat0, res0, clustrules)

% divide cluster into subregions.
stat0    = get_regions_kh(stat0, res0);

%%
if ops.splitROIs
    [stat, res] = get_validregions(stat0,res0, clustrules);
else
   stat = stat0;
   res = res0;
   
   for j = 1:length(stat)
      stat(j).igood = 1; 
   end
end


if ~exist(ops.ResultsSavePath, 'dir')
    mkdir(ops.ResultsSavePath)
end

Nk = ops.Nk;

% Nk      = numel(unique(res.iclust));

ops.ROISaveFile = sprintf('%s/F_%s_%s_plane%d_Ch%d_Nk%d.mat', ops.ResultsSavePath, ...
    ops.mouse_name, ops.date, ops.PlaneID,ops.ChannelID, Nk);
save(ops.ROISaveFile,  'ops', 'res', 'stat', 'stat0', 'res0', 'clustrules')

