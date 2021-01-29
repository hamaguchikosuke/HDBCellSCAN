function [Overlap,IA,IB] = ROIOverlapMatrix_001(ROIA,ROIB)
%% calculate overlaps of all the ROIs in the two images
% Note: the images need to be registered (aligned) so that 
% same ROI on different sessions are located close to each together.
% 
% Usage: [Overlap,IA,IB] = ROIOverlapMatrix_001(ROIA,ROIB)
% 
% ----- Input=-----
% ROIA, ROIB
% ROIA and ROIB is a matrix of integers where integers prepresents the ROI
% ID of each pixel. This program assumes that a single pixel belongs to a single ROI. 
% 
% Output: Overlap
% 
% Overlap = 2* #pixels in the intersection of ROI(1) and ROI(2)
%            -----------------------------------
%             #pixels in ROI(1) +  ROI(2)
% 
% IA: sorted labels in ROIA
% IB: sorted labels in ROIB.
%   
% Overlap(i,j) means overlap between IA(i) in ROIA and IB(j) in ROIB.
% 
% Overlap is larger than or equal to Jaccard distance which is defined as 
% Jaccard =   #pixels in the intersection of ROI(1) and ROI(2) 
%            -----------------------------------
%             #pixels of union( ROI(1) +  ROI(2))
% 
% 
% by KH 20170818
%%

% ROIA = MultiIMG(1).ROIs_numbered_moved;
% ROIB = MultiIMG(2).ROIs_numbered_moved;

% construct sparse matrix that represents each ROI.
ROIA_ID = sort(unique(ROIA));
ROIA_ID(ROIA_ID==0)=[];
ROIB_ID = sort(unique(ROIB));
ROIB_ID(ROIB_ID==0)=[];

IA= ROIA_ID;
IB= ROIB_ID;

ROIA_matrix = spalloc(prod(size(ROIA)),length(unique(ROIA))-1,nnz(ROIA));
for ii=1:length(ROIA_ID)
ROIA_matrix(find(ROIA==IA(ii)),ii)=1;
end

ROIB_matrix = spalloc(prod(size(ROIB)),length(unique(ROIB))-1,nnz(ROIB));
for ii=1:length(ROIB_ID)
ROIB_matrix(find(ROIB==IB(ii)),ii)=1;
end

Overlap = ROIA_matrix'*ROIB_matrix;
Npix = full(bsxfun(@plus,sum(ROIA_matrix,1)',sum(ROIB_matrix,1)));
Overlap = (2*Overlap./Npix);

