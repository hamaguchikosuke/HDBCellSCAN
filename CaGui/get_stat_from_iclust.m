function stat = get_stat_from_iclust(iclust,M)
% calculate stats from iclust and M (mean image) 
% stat = get_stat_from_iclust(iclust,M);
% --- input ---
% iclust: H x W interger matrix, each integer indicates the ROI index.
% M:     image of which pixel values are max correlation to neighboring pixel 
% 
% This function cannot be used when ROIs are overlapping, because iclust contains which pixel belongs to which ROI, 
% all or none. 
% 
% by KH 20170928

Ly = size(iclust,1);
Lx = size(iclust,2);

xs = repmat(1:Lx, Ly, 1);
ys = repmat((1:Ly)', 1, Lx);

xlx         = repmat(-ceil(Lx/2):1:ceil(Lx/2), 2*ceil(Lx/2)+1, 1);
rgrid       = sqrt(xlx.^2 + xlx'.^2);
rgridsort   = sort(rgrid(:), 'ascend');

Nk = numel(unique(iclust));
unique_clust = unique(iclust);
unique_clust(unique_clust==0)=[]; % in case zero is included as label, that is background.
for k = 1:max(unique_clust)%unique_clust(:)'
    ipix = find(iclust==k);
    if ~isempty(ipix)
        
        tmp = zeros(Ly,Lx);
        tmp(ipix)=1;
        
        x0 = xs(ipix); y0 = ys(ipix);
        if length(ipix)==0
            fprintf('cluster %d has length(ipix)=0\n');
        end
        edge_tmp = edge(tmp);
        
        sttmp = regionprops(tmp,'Eccentricity','Solidity','Perimeter');
        rs = ((x0 - median(x0)).^2 + (y0 - median(y0)).^2).^.5;    % spatial variance
        stat(k).mrs     = median(rs);                              % median value of spatial variance
        stat(k).npix    = numel(ipix);
        stat(k).mrs0    = median(rgridsort(1:stat(k).npix));
        stat(k).Compactness  = stat(k).mrs/stat(k).mrs0;
        stat(k).med     = [median(y0) median(x0)];
        stat(k).ipix    = ipix;                % index of labeled region
        stat(k).ipix_edge = find(edge_tmp);    % index of boundary pixels
        stat(k).lambda  = M(ipix);
        stat(k).V       = sum(M(ipix)); % added 20171102 to be consistent with get_regions.m
        stat(k).Solidity = sttmp.Solidity;   % npix/ConvexArea.
        stat(k).Eccentricity = sttmp.Eccentricity; % 0 is circle, 1 is line.
        stat(k).Perimeter = sttmp.Perimeter; % 0 is circle, 1 is line.
        stat(k).ipix_neuropil = []; % neuropil will be computed in get_signals_and_neuropil_kh.
        
    else
        stat(k).mrs     = NaN;                              % median value of spatial variance
        stat(k).npix    = 0;
        stat(k).mrs0    = NaN;
        stat(k).Compactness  = NaN;
        stat(k).med     = NaN;
        stat(k).ipix    = 0;                % index of labeled region
        stat(k).ipix_edge = NaN;    % index of boundary pixels
        stat(k).lambda  = NaN;
        stat(k).V       = NaN; % added 20171102 to be consistent with get_regions.m
        stat(k).Solidity = NaN;   % npix/ConvexArea.
        stat(k).Eccentricity = NaN; % 0 is circle, 1 is line.
        stat(k).Perimeter = NaN; % 0 is circle, 1 is line.
        stat(k).ipix_neuropil = []; % neuropil will be computed in get_signals_and_neuropil_kh.
    end
    
end