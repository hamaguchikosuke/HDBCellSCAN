function [L,probabilities,M,Sneuropil] = HDBCellScan_20171102(U,Sv,min_cluster_size)
% Usage: labels = HDBCellScan_20171102(U,Sv,min_cluster_size)
% === input ===
% U, Sv; eigenimage and power calculated from U*Sv*V = F, or [U,S,V]=svd(F); where F is Ca movie.
% min_cluster_size: about the minimum area of a cell. 
% 
% === output ===
% L: W x H matrix of labels (integers).
% 
% Usage2): 
% [L,prob,outlier] = HDBCellScan_20171102(U,Sv,min_cluster_size)
% 
% 
% ex) 
% % load('C:\home\ImagingExperiments\SimulatedData\SimulatedOne\20171020\1\SVDroi_SimulatedOne_20171020_plane1_ch1.mat');
% % load C:\home\ImagingExperiments\B6J647\SVDroi_B6J647_CREBRC107_20170902_d2_plane1_ch1.mat;
% load('C:\home\ImagingExperiments\B6J647\SVDroi_B6J647_CREBRC107_20170904_d6_plane1_ch1.mat','U','Sv');
% min_cluster_size = 70;
% 
% [L,probabilities]  = HDBCellScan_20171102(U,Sv,min_cluster_size);
% 
% [Ly,Lx]=size(L);
% L(L==0)=length(unique(L))+1;
% 
% % rng('default');
% rand_hue  = rand(1,length(unique(L)));
% rng('shuffle');
% rand_hue= [rand_hue,0];
% Hue = rand_hue(L);
% Sat = reshape(probabilities,Ly,Lx);
% Val = Sat;
% IMG=hsv2rgb(cat(3,Hue,Sat,Val));
% imagesc(IMG);
% 
% Usage3)
% [L,probabilities,M,S] = HDBCellScan_20171102(U,Sv,min_cluster_size);
% M: cos value of inner products (correlation level)  
% S: neuropil
% 
% 
% by Kosuke Hamaguchi 20171102

% Based on C:\home\Platex\2017Oct\test_hdbscan_CaImaging_20171028.m
% addpath('C:\home\matlab_svn\MatPy');

[Ly,Lx,nSVD]=size(U);
 
     
% myfigure('Obtained EigenImage');clf;
% NN=3;
% for ii=1:NN^2
%     mysubplot(NN,NN,ii); imagesc(U(:,:,ii));
%     box off;
%     axis off;
% end 

npix = Ly*Lx;
U = reshape(U,[],nSVD);
% U = bsxfun(@times, U, Sv'.^.5);
% U=U.*sqrt(Sv'); % > 2016b
U=U.*Sv';
% U=normc(U')';
%% in case to subtract neuropil (smooth components)

do_show_images = 0;
do_subtract_neuropil = 1; % 
% well, I need to do something to subtract neuropil for deeper imaging, or
% higher expression brain. 

if (do_subtract_neuropil)
    TileFactor = 1; % this option can be overwritten by the user
    cell_diameter = 20;
    nTiles = ceil(TileFactor * (Ly+Lx)/2 / (10 * cell_diameter)); % neuropil is modelled as nTiles by nTiles
    
    xc = linspace(1, Lx, nTiles);
    yc = linspace(1, Ly, nTiles);
    yc = yc';
    xs = 1:Lx;
    ys = 1:Ly;
    
    sigx = 4*(Lx - 1)/nTiles;
    sigy = 4*(Ly - 1)/nTiles;
    
    S = zeros(Ly, Lx, nTiles, nTiles, 'single');
    for kx = 1:nTiles
        for ky = 1:nTiles
            cosx = 1+cos(2*pi*(xs - xc(kx))/sigx);
            cosy = 1+cos(2*pi*(ys - yc(ky))/sigy);
            cosx(abs(xs-xc(kx))>sigx/2) = 0;
            cosy(abs(ys-yc(ky))>sigy/2) = 0;
            
            S(:, :,ky, kx) = cosy' * cosx;
        end
    end
    S = reshape(S, [], nTiles^2);
    S = normc(S);
    
    nBasis = size(S,2) ;
    PixL = ones(1, Lx * Ly)';
    
    Uneu = U';
    
    Sm = bsxfun(@times, S, PixL);
    StS = Sm' * Sm;
    StU = Sm' * Uneu';
    Lam = (StS + 1e-4 * eye(nBasis)) \ StU;
    
    % recompute neuropil pixel contribution
    neuropil = Lam' * S';
    PixL = mean(bsxfun(@times, neuropil, Uneu), 1);
    PixL = bsxfun(@rdivide, PixL, mean(neuropil.^2,1)); % normalize variance to 1
    PixL = max(0, PixL);
    neuropil = bsxfun(@times, neuropil, PixL);
    U = U - neuropil'; %what's left over for cell model
        
%     U = reshape(U,[],nSVD);
    fig_title = 'U-neuropil';
    Sneuropil = S;
else
    Sneuropil = [];
    fig_title = 'U-only';
end

if (do_show_images)
    figH=myfigure(fig_title);clf;
    NN=3;
    for ii=1:NN^2
        Utmp = reshape(U(:,ii),Ly,Lx);
        mysubplot(NN,NN,ii); imagesc(Utmp);
        box off;
        axis off;
    end
    drawnow;
end
    
%% construct neighbor indices
neigh = get_neighbor_index(Ly,Lx,'3x3','circular'); % I know circular boundary is not correct, but to avoid disconnected pixels at the corner, I used circular. 
% neigh = get_neighbor_index(Ly,Lx,'5x5','circular'); % I know circular boundary is not correct, but to avoid disconnected pixels at the corner, I used circular. 
%% check whether neigh vector contains neighboring index
do_check_neigh =0;

if (do_check_neigh)
    clf;
    % neigh(neigh==0)=length(neigh)+1;
    for ii=1:Ly*Lx
        A=zeros(Ly,Lx);
        ind = neigh(ii,find(neigh(ii,:)));
        A(ind)=1;
        imagesc(A); title(num2str(ii));
        pause(0.5);
    end
end
%% For each neighbor, calculate inner product of U.*(Sv.^0.5).
% To do this, construct sparse matrix which non-zero values are inner product of 
% feature vector U, but only between the 8-neighbors.
fprintf('Constructing sparse correlation matrix of nearest neighbors from svd eigenvectors...\n');
ind_i = neigh';
ind_j = repmat([1:Ly*Lx],size(neigh,2),1); % ind_j is column number.
ind_j=ind_j(:); ind_i=ind_i(:);
ind_j(ind_i==0)=[]; 
ind_i(ind_i==0)=[];

nbatch = Ly*Lx; 
nind = length(ind_i);
start_ind = 1:nbatch:nind;
end_ind   = start_ind+nbatch-1; end_ind(end)=nind;
s = zeros(1,length(ind_i));
cnt=0;
nU =normr(U);

for ii=1:length(start_ind)
    index = start_ind(ii):end_ind(ii);
    s(cnt+[1:length(index)])=sum(nU(ind_j(index),:).*nU(ind_i(index),:),2)';
    cnt = cnt+length(index);
end

% construct sparse matrix which indicates the 8-neighbors connection.
s1 = s;         
s = max(s)-s+eps; % make it distance
mins = min(s);
maxs= max(s);
s = (s/(maxs-mins)); % normalize to 0 to 1

S = sparse(ind_i,ind_j,s,Ly*Lx,Ly*Lx);
M = sparse(ind_i,ind_j,s1,Ly*Lx,Ly*Lx);
M = full(max(M,[],1)); % most correlated pixels (max values)

%% Finally, test with hdbscan: well it did not work. Finally it worked. 
pyS = matsparse_2_pysparse(S);
 tic
 fprintf('HDBCellSCAN...(it will take less than a minutes (if Python 3.5>))');
 
%  import sys
%     sys.path.append('C:\\Users\\hamag\\.conda\\envs\\khtest\\Lib\\site-packages')
%     import hdbscan
    PA = pyargs('min_cluster_size',uint16(min_cluster_size),...
                'min_samples', uint16(2),...
                'gen_min_span_tree',true, ...
                'metric','precomputed',...
                'cluster_selection_method','leaf',...
                'alpha',1);
%              'cluster_selection_method','leaf',...
    clusterer = py.hdbscan.HDBSCAN(PA);
    clusterer=clusterer.fit(pyS);
    fprintf('Done\n')
    toc
% clusterer=py.my_hdbscan.hdbSfit(pyS); % this is the function to do hdbscan 
 labels=nparray2mat(clusterer.labels_)+2; % In pythonm -1 is noise, cluster starts from 0.
%  Now in matlab, 1 is noise, cluster starts from noise cluster 1.
 probabilities=nparray2mat(clusterer.probabilities_);
 
 %%
 L=reshape(labels,Ly,Lx);
 % 20171115: changed. iclust==0 means background. 
%  L(L==0)=length(unique(labels)); % last label is for background.
 L_tmp = L;
  L_tmp(L_tmp==0)=length(unique(labels)); 
 M = reshape(M,Ly,Lx);
 
 if (do_show_images)
     figure(figH);clf;
     % rng('default');
     rand_hue  = rand(1,length(unique(labels)));
     rng('shuffle');
%      rand_hue= [0 rand_hue,0];
     Hue = rand_hue(L_tmp);
     Sat = reshape(probabilities,Ly,Lx);
     Val = Sat;
     IMG=hsv2rgb(cat(3,Hue,Sat,Val));
     imagesc(IMG);
     drawnow;
 end
%% in case some cells are merged, I need to develop a method to separate them.
 % obtain left-top and right-bottom corner.
%  x=ceil(ginput(2));
%  J= x(1,2):x(2,2);
%  I= x(1,1):x(2,1);
%  index=reshape(1:Ly*Lx,Ly,Lx);
%  index=index(J,I);
% %  Utmp = U(index(:),:);
%  Stmp =S(index(:),index(:));
%  
%  pyS = matsparse_2_pysparse(Stmp);
%  tic
%     PA = pyargs('min_cluster_size',uint32(70),...
%                 'min_samples', uint32(2),...
%                 'gen_min_span_tree',true, ...
%                 'metric','precomputed',...
%                 'cluster_selection_method','eom',...
%                 'alpha',1);
%     clusterer2 = py.hdbscan.HDBSCAN(PA);
%     clusterer2=clusterer2.fit(pyS);
%    
%  labels2=nparray2mat(clusterer2.labels_)+1; % python index starts from 0.
%  probabilities2=nparray2mat(clusterer2.probabilities_);
%  
%  myfigure('PCA of Utmp');clf;
%  L2=reshape(labels2,length(J),length(I));
% L2(L2==0)=length(unique(labels2))+1;
% rand_hue  = rand(1,length(unique(labels2)));
% rand_hue= [rand_hue,0];
% Hue = rand_hue(L2);
% Sat = reshape(probabilities2,length(J),length(I));
% Val = Sat;
% IMG=hsv2rgb(cat(3,Hue,Sat,Val));
% imagesc(IMG);
% 
% myfigure('dendrogram of Utmp')
% slt=nparray2mat(clusterer2.single_linkage_tree_.to_numpy);
% dendrogram([slt(:,1:2)+1,slt(:,3)]); % first two columns are index starting from 0, 3rd is the distance.
% 
%% well, dendrogram did not reveal the diffrence. Let's check the distribution of pixel.
% 
% % now, pickup the two points which you want to divide.
% x = ceil(ginput(2));
% 
% %
% % if iclust1~=iclust2
% %     error('Please select the same ROI to split');
% % end
% index=reshape(1:Ly*Lx,Ly,Lx);
% indY=[];
% indX=[];
% indAll = [];
% seed_ind = zeros(size(x,1),1);
% seed_ind_in_indAll = zeros(size(x,1),1);
% iclust=zeros(size(x,1),1);
% Val = zeros(Ly,Lx);
% 
% for ii=1:size(x,1)
%     tmp=L(x(ii,2),x(ii,1));
%     ind=find(L==tmp);
%     fprintf('iclust%d=%d\n ',ii,tmp);
%     seed_ind(ii)=index(x(ii,2),x(ii,1));
%     seed_ind_in_indAll(ii)=find(ind==seed_ind(ii));
%     if any(iclust==tmp)
%         % detected same cluster
%         fprintf('Detected the same cluster, skip!\n');
%         iclust(ii)=[];
%         continue;
%     else
%         iclust(ii) = tmp;
%     end
%     
%     % 
%     
%     % get the square region 
%     [index_y,index_x]=ind2sub([Ly,Lx],ind);
%     indY=cat(1,indY,index_y);
%     indX=cat(1,indX,index_x);
%     indAll{ii} = ind;
%     Val = Val | L==iclust(ii);
%     
%    
% end
% J=min(indY):max(indY);
% I=min(indX):max(indX);
% myfigure('Utmp plot');clf;
% subplot(2,2,1);
% 
% Hue = rand_hue(L);
% Sat = reshape(probabilities,Ly,Lx);
% IMG=hsv2rgb(cat(3,Hue,Sat,Val));
% imagesc(IMG(J,I,:));
% 
% 
% % plot the principle components.
% all_index = cat(1,indAll{:});
% Utmp = U(all_index,:);
% Utmp2 = cat(2,zscore(indY),zscore(indX),zscore(Utmp,0,1)/4); % 
% [coef,score] = pca(Utmp2);
% subplot(2,2,2);cla;
% plotH=[];
% ind=0;
% Col = {'r','g','b'}; 
% for ii=1:length(indAll)
%     ind=ind+[1:length(indAll{ii})];
%     plotH(ii)=plot3(score(ind,1),score(ind,2),score(ind,3),[Col{ii},'.']);hold on;
%     indtmp=ind(seed_ind_in_indAll(ii));
% %     [score(indtmp,1), score(indtmp,2), score(indtmp,3)];
%     plot3(score(indtmp,1), score(indtmp,2), score(indtmp,3),[Col{ii},'o']);
%         ind = ind(end); 
% %       pause
% end
% 
% grid on;
% %  next, plot the result of k-means
% 
% Nlabel = 2;
% 
%  idx = kmeans(Utmp2,Nlabel,'Start',Utmp2(seed_ind_in_indAll,:));
%  subplot(2,2,3);cla;
%  for ii=1:Nlabel
%      ind=find(idx==ii);
%      plotH(ii)=plot3(score(ind,1),score(ind,2),score(ind,3),[Col{ii},'.']);hold on;
%  end
% grid on;
% 
% % recover the divided image 
% 
% subplot(2,2,4);cla;
% 
% Hue = rand_hue(L);
% Sat = reshape(probabilities,Ly,Lx);
% 
% Col = [0.9,0.4,0.8];
% % Val = L~=0;
% 
% for ii=1:Nlabel
%     original_ind = all_index(idx==ii);
%     Hue(original_ind)=Col(ii);
% %     Sat(original_ind)=1;
% end
% 
% IMG=hsv2rgb(cat(3,Hue,Sat,Val));
% image(IMG(J,I,:));
% 
% %% Find the cell bodies by morphological opening 
% 
% L0=L; L0(L==max(L(:)))=0; % make noise as zero.
% [B,CoM,shrinkL0] = labeledIMG_bwtraceboundary(L0,3);
% 
% % se=strel('disk',1,4);
% % mask=L; mask(L==max(L(:)))=0; % make noise as zero.
% % mask = mask~=0;
% % Nrepeat = 5;
% % for nn=1:Nrepeat,    mask = imerode(mask,se,'same'); end
% % for nn=1:Nrepeat,    mask = imdilate(mask,se,'same'); end
% 
% shrinkL = shrinkL0;
% shrinkL(shrinkL==0)=max(L(L(:)));
% myfigure('cell bodies');clf;
% Hue = rand_hue(shrinkL);
% Sat = reshape(probabilities,Ly,Lx);
% Val = shrinkL0~=0;
% 
% IMG2=hsv2rgb(cat(3,Hue,Sat,Val));
% image(IMG2);hold on;
% 
% for ii=1:length(B)
%   xy=B{ii};
%   if ~isempty(xy)
%       plot(xy(:,2),xy(:,1),'Color',[1 1 1]); hold on;
%   end
% end
% 
% %% or we can try filter 
% figh=myfigure('ROI check');cla;
% ax(1)=subplot(2,2,1);
% L0=L; L0(L==max(L(:)))=0; % make noise as zero.
% L0=double(L0~=0);
%  se=strel('disk',3,4);
%  L0= imopen(L0,se);
% imagesc(L0);
% % colormap gray;
% 
% 
% celldiameter = 11;
% cell_membrane = 2;
% cell_edge  = 2;
% cellR = celldiameter/2;
% 
% % L0=zeros(101,101);
% % L0(50,50)=1;
% 
% L0_cellbody = imgaussfilt(L0,cellR*1.2);
% L0_boundary = imgaussfilt(L0,cellR*1.9);
% L0_cell = L0_cellbody-L0_boundary;
% subplot(2,2,2);cla;
% imagesc(L0_cell);
% 
% %
% subplot(2,2,3);cla;
% localMin = ordfilt2(L0_cell,1,ones(celldiameter));
% L0_cell = L0_cell-ordfilt2(localMin,celldiameter^2,ones(celldiameter));
% imagesc(L0_cell);
% 
% localMax = ordfilt2(L0_cell,celldiameter^2,ones(celldiameter));
% imagesc(localMax);
% L0peaks = (L0_cell>localMax-eps) &( L0>0);
% imagesc(L0+2*L0peaks);
% %%
% hsize = celldiameter+cell_edge+cell_membrane;
% h_boundary = zeros(hsize,hsize); 
% h_boundary_ind = 1:hsize;
% h_boundary(h_boundary_ind,h_boundary_ind)=fspecial('disk',(length(h_boundary_ind)-1)/2);
% 
% h_nuclei       = zeros(hsize,hsize);
% h_nuclei_ind = (1+cell_membrane+cell_edge):(hsize-cell_membrane-cell_edge);
% h_nuclei(h_nuclei_ind,h_nuclei_ind)=fspecial('disk',(length(h_nuclei_ind)-1)/2);
% 
% h_mem = zeros(hsize,hsize); 
% h_mem_ind = (1+cell_membrane):(hsize-cell_membrane);
% h_mem(h_mem_ind,h_mem_ind)=fspecial('disk',(length(h_mem_ind)-1)/2);
% 
% subplot(2,2,2);
% h = 2*h_mem-h_nuclei-h_boundary;
% imagesc(h);
% 
% subplot(2,2,3);
% conL0= conv2(L0,h,'same');
% ax(2)=subplot(2,2,2);
% imagesc(conL0);
% 
% %
% 
% nimg_rgb = repmat(mat2gray(nimg),[1,1,3]);
% cell_center = zeros(size(nimg_rgb));
% cell_center(:,:,1)=conL0>0.8;
% cell_center(:,:,2)=conL0<-0.9;
% ax(3)=subplot(2,2,3);
% imagesc(tanh(nimg_rgb+cell_center));
% 
% 
% % deconv_cell_center=zeros(size(nimg_rgb));
% % deconv_cell_center(:,:,1)=imdeblur(double(cell_center(:,:,1)),fspecial('gaussian',[5 5],1));
% % ax(4)=subplot(2,2,4);
% % imagesc((1+tanh(nimg_rgb+deconv_cell_center))/2);
%  
% linkaxes(ax,'xy');

