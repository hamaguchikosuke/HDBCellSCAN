function inspect_U(U,ind,varargin)
% inspect_U(U,1:10)
% inspect_U(U,1:10,[Ly,Lx]);

if ndims(U)==2
    if nargin>=3
        Lyx=varargin{1};
        Ly=Lyx(1);
        Lx=Lyx(2);
        U=reshape(U,Ly,Lx,[]);
    else
        error;
    end
end

myfigure(mfilename);clf;
for ii=1:length(ind)
    imagesc(U(:,:,ind(ii)));
    title(num2str(ind(ii)))
    pause(1);
end
    