function stat = get_stat2(res)

Ly = res.Ly;
Lx = res.Lx;

xs = repmat(1:Lx, Ly, 1);
ys = repmat((1:Ly)', 1, Lx);

xlx         = repmat(-ceil(Lx/2):1:ceil(Lx/2), 2*ceil(Lx/2)+1, 1);
rgrid       = sqrt(xlx.^2 + xlx'.^2);
rgridsort   = sort(rgrid(:), 'ascend');

Nk = numel(unique(res.iclust));
for k = 1:Nk
%     ipix = find(res.iclust==k);    
%     x0 = xs(ipix); y0 = ys(ipix);
    [y0,x0] = find(res.iclust==k);    
    
    rs = ((x0 - median(x0)).^2 + (y0 - median(y0)).^2).^.5;    % spatial variance    
    stat(k).mrs     = median(rs);                              % median value of spatial variance
    stat(k).npix    = numel(x0);
    stat(k).mrs0    = median(rgridsort(1:stat(k).npix));
    stat(k).med     = [median(y0) median(x0)];
    stat(k).ipix    = ipix;
    stat(k).lambda  = res.M(ipix);
    


%     stat(k).iscell = stat(k).mrs/stat(k).mrs0<1.3 & ...
%         stat(k).npix>50 & stat(k).npix<300;    
end