function [tiffname,tiffpaths]=get_tiff_location_from_ops_20171102(ops,varargin)

% this ops needs to be build by 
% ops = build_ops3(ops) 
% otherwise, to it here.
if isfield(ops,'ResultsSavePath')
    % OK
else
    ops = build_ops3(ops);
end

if nargin>=2
    PlaneChString = varargin{1};
else
    PlaneChString = ops.PlaneChString;
end

tiffname = sprintf('%s_%s_%s_2P_%s_%s.tif', ops.date, '*', ...
    ops.mouse_name, 'plane*_ch*', '*');

for kk=1:length(ops.SubDirs)
tiffpaths{kk} =  fullfile(ops.RegFileTiffLocation, ops.mouse_name, ops.date, ...
            ops.SubDirs{kk}, PlaneChString);
end