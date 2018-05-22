 function [ops, stat, res] = fast_spect_clustering_KH(ops, U, Sv)
% k-means clustering of pixels in feature space U. 
% by Kosuke Hamaguchi, largely based on Pichatariau's fast clustering
% (which uses inner products of centroid vector and pixel's feature).
% Here, I used distance in feature space U. 
% 20170308

U =  reshape(U, [], size(U,ndims(U)));
iplane = ops.iplane;

for i = 1:size(U,2)
  U(:,i) = U(:,i)  * Sv(i).^.5;
end
U = U';
[nSVD, Npix] = size(U);
%% clustering options
Ly = numel(ops.yrange);
Lx = numel(ops.xrange);

ops.Nk0 = ceil(sqrt(ops.Nk0)).^2;

NCluster      = ops.Nk0;
niter   = ops.niterclustering;



xs = repmat(1:Lx, Ly, 1);
ys = repmat((1:Ly)', 1, Lx);



randx = rand(1, NCluster) * Lx;
randy = rand(1, NCluster) * Ly;

dx = repmat(xs(:), 1, NCluster) - repmat(randx, numel(xs(:)), 1);
dy = repmat(ys(:), 1, NCluster) - repmat(randy, numel(ys(:)), 1);

dxy = dx.^2 + dy.^2;

%%
[~, iclust] = min(dxy, [], 2);

clear dx dy

if ops.ShowCellMap
    myfigure('CellMap');
%     figure( 'Units', 'pixels', 'position', [100 100 900 900])
    colormap('hsv')
    axes('position', [.05 .05 .925 .925])
    set(gcf, 'Color', 'w')
end

r   = rand(1, NCluster);
Sat = ones(Ly, Lx);

err = zeros(niter,1);
ops.meanV = gather(sum(Sv)/(Ly*Lx));

NCluster = ops.Nk;
NClusterIter = round(linspace(ops.Nk0, ops.Nk, niter-2));
NClusterIter(end+1:(niter+1)) = ops.Nk;

M = ones(Npix,1, 'single');
%%
tic
ison = true(NCluster,1);
NoiseLevel = 0.001;
d = zeros(NCluster,Npix);
NInternalRepeats = 2;
c = zeros(nSVD, NCluster, 'single');
  
for k = 1:niter  
%%     [NClusterIter(k) numel(unique(iclust)) sum(ison)]


for cc=1:NInternalRepeats
    % update centroid position
    for kk=1:NCluster
        ind = iclust==kk;
        if isempty(ind)
            c(:,kk)=nan(nSVD,1);
        else
            c(:,kk)=mean(U(:,ind),2);
        end
    end
    
    
    % calculate distance to assign cluster id to each pixel
    for kk=1:NCluster
        tmp=bsxfun(@minus,U,c(:,kk));
        d(kk,:)=sqrt(sum(tmp.^2,1));
        %         d(kk,:)=sqrt(sum(tmp.^2,1))+NoiseLevel*rand(1,length(tmp)); % in some case, better to add some noise for decision.
    end
end
%
    
    % find the centroid that is closest to each pixel.
    [M, iclust] = min(abs(d),[],1);    
    d(iclust + (0:NCluster:numel(d)-1)) = inf;
    [M2, iclust2] = min(abs(d),[],1);
    
    dM = M - M2;
    indrem = NClusterIter(k) - NClusterIter(k+1);
    dMk = Inf*ones(NCluster,1);
    icl = cell(NCluster,1);
    
    for ii = 1:NCluster
        icl{ii} = iclust==ii;
        if ~isempty(icl{ii})
            dMk(ii) = sum(dM(icl{ii}));
        end
    end
    
    vlk = zeros(NCluster, 1);
    vlk(~ison) = Inf;
    while indrem>0
        % merge those pixels which centroid is too close to each other.
       [Xmin, imin]          = min(dMk + vlk);
       if isinf(Xmin)
          NClusterIter(k+1) = sum(ison);
          break; 
       end
       newi               = iclust2(icl{imin});
       iclust(icl{imin})  = newi;
       M(icl{imin})       = M2(icl{imin});
       dMk(unique(newi))  = Inf;
       dMk(imin)          = Inf;
       dMk(unique(iclust(iclust2==imin))) ...
           = Inf;
       
       ison(imin) = 0;
       indrem             = indrem - 1;
    end
    
    
    M       = M';
    iclust  = iclust';
 %%   
    err(k) = sum(M(:));
    if rem(k,2)==1 && ops.ShowCellMap
        %%
        lam = 1/M;
        for i = 1:NCluster
            ix = find(iclust==i);
            nT0 = numel(ix);
            if nT0>0
                vM = lam(ix);
                vM = vM/sum(vM.^2)^.5;
                lam(ix) = vM;
            end
        end
        V = max(0, min(.5 * reshape(lam, Ly, Lx)/mean(lam(:)), 1));
        H = reshape(r(iclust), Ly, Lx);
        rgb_image = hsv2rgb(cat(3, H, Sat, V));
        imagesc(rgb_image)
        title(sprintf('Iter=%d',k));
        axis off
        drawnow
        fprintf('explained variance is %2.6f time %2.2f \n', err(k), toc)
    end
   
    if  sum(ison)<NCluster
        fprintf('Converged.\n')
        break; 
    end
end

lam = 1/M;
for i = 1:NCluster
    ix = find(iclust==i);
    
    nT0 = numel(ix);
    if nT0>0
        vM = lam(ix);
        lam(ix) = vM/sum(vM.^2)^.5;
    end
end

%%
newindx = cumsum(ison);
iclust  = newindx(iclust);
NCluster      = numel(unique(iclust));
%%
clear res

res.iclust  = iclust;
res.M       = M;
res.lambda  = lam;

%%
res.Ly  = Ly;
res.Lx  = Lx;
stat    = get_stat(res);
%%
if ~exist(ops.ResultsSavePath, 'dir')
    mkdir(ops.ResultsSavePath)
end
save(sprintf('%s/F_%s_%s_plane%d_Nk%d.mat', ops.ResultsSavePath, ...
    ops.mouse_name, ops.date, iplane, NCluster),  'ops', 'res', 'stat')

%%
% sk = skewness(F,[],2);
% [~, isk] = sort(sk, 'descend');
% clf
% for i = 1:20
%    plot(5*i + zscore(F(isk(i), :)))
%    hold all
% end
% axis tight
