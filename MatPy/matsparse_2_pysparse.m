function pyS = matsparse_2_pysparse(X)
% Convert matlab sparse matrix to python scipy sparse format.
% Usage: pyS = matsparse_2_pysparse(X)
% X is the sparse matrix.
% 
% --- here is the basic python code to make a sparse matrix ---
% import scipy.sparse 
% from scipy.sparse import csr_matrix
% from scipy import array
% 
% # create a sparse matrix
% row = array([0,0,1,2,2,2])
% col = array([0,2,2,0,1,2])
% data = array([1,2,3,4,5,6])
% 
% mat = csr_matrix( (data,(row,col)), shape=(3,4) )
% 
% # get the data, shape and indices
% (m,n) = mat.shape
% s = mat.data
% i = mat.tocoo().row
% j = mat.indices
% 
% # display the matrix
% print mat
% --- and here is a code to convert a matlab sparse matrix to scipy sparse matrix 
% m = 3;
% n = 4;
% s = [1, 2, 3, 4, 5, 6];
% %Index from 1 in Matlab.
% i = [0, 0, 1, 2, 2, 2] + 1;
% j = [0, 2, 2, 0, 1, 2] + 1;
% 
% S = sparse(i, j, s, m, n, m*n);
% % in Python, () is tuple, []
%  
% inputargs = pyargs('shape',py.tuple({m,n}));
% s = py.scipy.array(s);
% i = py.scipy.array(uint32(i-1));
% j = py.scipy.array(uint32(j-1));
% input1 = py.tuple({s,py.tuple({i,j})}); % equivalent to (s, (i,j)) in python
% pyS=py.scipy.sparse.csr_matrix(input1,pyargs('shape',py.tuple({m,n})));
% 
% % compare 
% full(S)
% nparray2mat( pyS.toarray)
%%          

if issparse(X)
    % OK
else
    error('Input X must be a sparse matrix');
end

[I,J,S]=find(X);
inputargs = pyargs('shape',py.tuple({py.int(size(X,1)),py.int(size(X,2))}));
S = py.scipy.array(S');
I = py.scipy.array(uint32(I'-1));
J = py.scipy.array(uint32(J'-1));
input1 = py.tuple({S,py.tuple({I,J})}); % equivalent to (S, (I,J)) in python
pyS=py.scipy.sparse.csr_matrix(input1,inputargs);
