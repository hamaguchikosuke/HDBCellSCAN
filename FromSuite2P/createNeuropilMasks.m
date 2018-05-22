function neuropMasks=createNeuropilMasks(cellFields,allField,xPU,yPU,ops)
%
% createNeuropilMasks compute the a ring-like shape around a cell mask to
% determine the neuropil contamination. The mask excludes any mask
% belonging to other cells.
%
% INPUTS:
% cellFields: binary [nCells x nY x nX] matrix containining cell masks
% allField: [nY x nX] matrix identify different cell masks with an integer
% index
% xPU: pixels per micrometers (horizontal axis)
% yPU: pixels per micrometers (vertical axis)
% ops.inNeurop: distance between cell mask and inner bondary of neuropil
% mask (in um), default is 3 um
% ops.outNeurop: outer radio of Neuropil mask (in um), default is 20 um
%
% OUTPUTS:
% allNeurop: [nY x nX] matrix identify different neuropil masks with an
% integer index
%
% 2015.06.11 Mario Dipoppa - Created
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3 || ~isfield(ops, 'inNeurop')
    inNeurop=3;
else
    inNeurop=ops.inNeurop;
end
if nargin<3 || ~isfield(ops, 'outNeurop')
    outNeurop=20;
else
    outNeurop=ops.outNeurop;
end

inRadius=ceil(inNeurop*sqrt(xPU*yPU));
outRadius=outNeurop*sqrt(xPU*yPU);

[hp,wp]=size(allField);


[xx_np, yy_np] = meshgrid((1:wp),(1:hp));

nCells=size(cellFields,1);

%compute the inner border of the ring-like shaped surrounding neuropils
se = strel('disk',inRadius);
expandedGeneralMask=sign(imdilate(allField,se));

neuropMasks=zeros(nCells,hp,wp);

%compute the mask of each surrounding neuropils
for ii=1:nCells
    stat = regionprops(squeeze(cellFields(ii,:,:)),'centroid');
    if isempty(stat)
        continue
    end
    centerCell=round(stat.Centroid);
    
    neuropCircle = sqrt((xx_np-centerCell(1)).^2+(yy_np-centerCell(2)).^2)<=outRadius;
    neuropMasks(ii,:,:)=(neuropCircle-expandedGeneralMask).*((neuropCircle-expandedGeneralMask)>0);
end

