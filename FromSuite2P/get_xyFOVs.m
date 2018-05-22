function [xFOVs, yFOVs] = get_xyFOVs(ops,varargin)

if nargin>=2
    RegShape = varargin{1};
else
    RegShape ='Square';
end

Ly = ops.Ly;
Lx = ops.Lx;

switch RegShape
    case 'Square'
        iR = ops.splitFOV;
        
        ny = floor(Ly/iR(1));
        nx = floor(Lx/iR(2));
        
        xFOVs = zeros(nx, iR(2), iR(1));
        yFOVs = zeros(nx, iR(2), iR(1));
        for i = 1:iR(1)
            for j = 1:iR(2)
                xFOVs(:,j,i) = [1:nx] + (j-1)*nx;
                yFOVs(:,j,i) = [1:ny] + (i-1)*ny;
            end
        end
        
        xFOVs = xFOVs(:,:);
        yFOVs = yFOVs(:,:);

    case 'GrinLens'
        % use center 2x2 from 4x4 split view. 
        iR = [4 4]; 
        
        ny = floor(Ly/iR(1));
        nx = floor(Lx/iR(2));
    
        xFOVs = nx+[1:2*nx]';
        yFOVs = ny+[1:2*ny]';
    otherwise
        error('Unknown option %s',RegShape);
end
