%  test for writing hdbscan code
 addpath C:\home\Platex\2017Oct\
 
 toolbox_path = 'C:\home\matlab_svn\matlab_bgl';
if exist(toolbox_path, 'dir')
	addpath(genpath(toolbox_path)) % add local path to the toolbox
else
	error('toolbox_path does not exist, please change toolbox_path');
end


 load clr-24-1.mat
%  spy(A)
 view(biograph(A,[],'ShowArrows','off','ShowWeights','on'));
 
 %% calculate minimum spanning tree  (requires matlab_bgl by David Gleich)
 tree = mst(A);
 view(biograph(tree,[],'ShowArrows','off','ShowWeights','on'));
 
%% now I want to convert it to dendrogram 
% linkage function requires pdist generated weight vector,
% which is ordered in (2,1), (3,1),..., (m,1), (3,2),(4,2),.... (index of lower triangle part of distance matrix)
Atmp = A;
Atmp(Atmp==0)=max(A(:))+10;
 triu_ind = find(tril(ones(size(Atmp)),-1));
 pdist_vecror = full(Atmp(triu_ind))';

 [dendro,N] = linkageD(pdist_vecror,'single');
clf;subplot(2,1,1);
dendrogram(dendro);
title('dendrogram from A');
 %% we should get the same results from minimum spanning tree when method = 'single'
%  because linkage analysis converts A to tree inside the code.
Atmp = tree;
Atmp(Atmp==0)=max(A(:))+10;
 triu_ind = find(tril(ones(size(Atmp)),-1));
 pdist_vecror = full(Atmp(triu_ind))';

 dendro = linkageD(pdist_vecror,'single');
 subplot(2,1,2); dendrogram(dendro);
title('dendrogram from minimum spanning tree of A')
 
 %%
 rmpath addpath C:\home\Platex\2017Oct\
 
 %% Now let's try using python hdbscan
 
 clear classes
 mod = py.importlib.import_module('my_hdbscan')
 py.importlib.reload(mod)
 %%
 
 Use_Sparse_Distance = 1;
 
 moons=py.sklearn.datasets.make_moons(pyargs('n_samples',50,'noise',0.05));
 % in Python, centers=[(-0.75,2.25), (1.0, 2.0)] is equal to 
 centers=py.list({py.tuple({-0.75, 2.25}), py.tuple({1.0, 2.0})});
 
 blobs=py.sklearn.datasets.make_blobs(pyargs('n_samples',int32(50),'centers',centers,'cluster_std',0.25));
 test_data = py.numpy.vstack(py.list({moons.cell{1}, blobs.cell{1}}));
 
  %  clusterer=py.my_hdbscan.fit3([]); % just get clusterer.
 clusterer=py.my_hdbscan.hdbfit(test_data); % this is the function to do hdbscan 
 labels=nparray2mat(clusterer.labels_)+1; % python index starts from 0.
 probabilities=nparray2mat(clusterer.probabilities_);
 
 
 test_data=nparray2mat(test_data);
 ColorPalette = [
     0         0    1.000;...       % b
     0     1.000        0;...       % g
     1.0000         0        0;...       % r
     0     1.000    1.000;...h
     1.0000         0    1.000;...a
     1.0000     1.000       	0;...
     0.7500    0.7500   0.2500;...
     0.7500    0.2500   0.7500];             %
 
 h = zeros(size(test_data,1),1);
 clf;
 for ii=1:size(test_data,1)
     if labels(ii)==0
         col = [0.7 0.7 0.7];
     else
         col=ColorPalette(labels(ii),:);
     end
     h(ii)=plot(test_data(ii,1),test_data(ii,2),'o','MarkerEdgeColor',col);hold on;
 end
 
 %% now try sparse matrix input

%  py_tree = matsparse_2_pysparse(tree);
 py_tree = matsparse_2_pysparse(A);
 
 clusterer=py.my_hdbscan.hdbfit(py_tree); % this is the function to do hdbscan 
 labels=nparray2mat(clusterer.labels_)+1; % python index starts from 0.
 probabilities=nparray2mat(clusterer.probabilities_);
  
%% well, it does not work. I should try more larger data set


T = zeros(Ly+2, Lx+2);
T(2:end-1, 2:end-1) = reshape(1:(Lx*Ly), Ly, Lx);
neigh = cat(3, T(1:end-2,2:end-1), T(3:end,2:end-1), T(2:end-1,1:end-2), T(2:end-1,3:end));
neigh = reshape(neigh, [], 4);
 
