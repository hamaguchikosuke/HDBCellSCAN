function index = get_neighbor_index(Ly,Lx,varargin)
% return index of neighbors in Lx x Lx matrix
% index = get_neighbor_index(Ly,Lx);
% 
% Usage: index = get_neighbor_index(Ly,Lx,varargin)
% === input ===  
% Ly: pixel number in y-axis; 
% Lx: pixel number in x-axis; 
% ex) 
% Ly =5;
% Lx =5;
% index = get_neighbor_index(Ly,Lx);
% 
% === output ===
% index: index of non-zero elements in distance matrix.
% 
% Usage2) index = get_neighbor_index(Ly,Lx,opt)
% opt: '3x3' is 8-nearest neighbors. % default
%    : '5x5' is 15-nearest neighbors. 
% 
% Ly =5;
% Lx =5;
% index = get_neighbor_index(Ly,Lx);
% ind_i = index';
% ind_j = repmat([1:Ly*Lx],size(index,2),1); % ind_j is column number.
% ind_j=ind_j(:); ind_j(ind_i==0)=[]; 
% ind_i=ind_i(:); ind_i(ind_i==0)=[];
% S = sparse(ind_i,ind_j,ones(1,length(ind_i)),Ly*Lx,Ly*Lx);
% spy(S);
% 
% Usage3) index = get_neighbor_index(Ly,Lx,opt,boundary);
% boundary : 'non-circular'; % default.
%          : 'circular';     % circular boundary condition.
% 
% by KH 20171020

if nargin<=2
    opt = '3x3';
else
    opt = varargin{1};
end

if nargin<=3
    boundary = 'non-circular';
else
    boundary = varargin{2};
end

switch opt
    case '3x3'
        nshift = 1;
    case '4x4'
        error('No 4x4 option!');
    case '5x5'
        nshift = 2;
    otherwise
        error('Unknown option %s',opt);
end
%% an easy way to construct neighbor indices
switch boundary
    case 'non-circular'
        T = zeros(Ly+2*nshift, Lx+2*nshift);
        T((1+nshift):end-nshift, (1+nshift):end-nshift) = reshape(1:(Lx*Ly), Ly, Lx);
        index_y = (1+nshift):(Ly+2*nshift-nshift);
        index_x = (1+nshift):(Lx+2*nshift-nshift);
    case 'circular'
        T = reshape(1:(Lx*Ly), Ly, Lx);
        index_y = 1:Ly;
        index_x = 1:Lx;  
    otherwise
        error('Unknown boundary condition %s!',boundary);
end


shifts_x = [nshift:-1:-nshift]; 
shifts_y = [nshift:-1:-nshift];
[shifts_X, shifts_Y]=meshgrid(shifts_x,shifts_y);
neigh = zeros(Ly,Lx,numel(shifts_X)-1);

cnt = 1;
for ii=1:numel(shifts_X)
    dx = shifts_X(ii);
    dy = shifts_Y(ii);
    if dx==0 && dy==0
        continue;
    else
%         [dy,dx]
        tmp = circshift(T,[dy,dx]);
        neigh(:,:,cnt)=tmp(index_y,index_x);
        cnt = cnt+1;
    end
end
index = reshape(neigh, [], size(neigh,3));

%% more straightforward way
% if nargin<3
%     index = 1:Ly*Lx;
% else
%     index = varargin{1};
% end
% index = index(:); % make it column vector;
% 
% mode = '8_neighbors';
% % npix = Ly*Lx;
% switch mode
%     case '8_neighbors'
%         shifts = [-Ly-1, -Ly, -Ly+1, -1, 1, Ly-1, Ly, Ly+1];
%     case '4_neighbors'
%         shifts = [-Ly, -1, 1,  Ly];
%     otherwise
%         error('Unknown mode %s',mode);
% end
% 
% index = bsxfun(@plus, index,shifts);
% index(index>Ly*Lx)=0;
% index(index<1)=0;