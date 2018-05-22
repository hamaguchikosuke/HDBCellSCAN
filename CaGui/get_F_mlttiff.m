function F = get_F_mlttiff(Ipix, FullTiffNames,varargin)
%  get mean fluorescent signals from multiple Tiffs 
%  Usage1) F = get_F_mlttiff(Ipix, FullTiffNames)
% 
%  --- Input1 ---
%  Ipix: 1xN or Nx1 cell array which contains indices of N- ROI. 
%  FullTiffNames: a cell of multiple tiff file names including all the path
% 
%  Usage2) F = get_F_mlttiff(Ipix, FullTiffNames,x_range,y_range)
% 
%  --- Input2 ---
%  x_range, y_range: Ipix is (often the case) indices in registered movie, 
%  which is smaller than the original Tiff file. To define the boundary of
%  the registered image, define x_range and y_range as 1 x W, 1xH vector.
%  If they are empty, x_range is assume to be 1:Width, and y_range = 1:Height of tiff file. 
% 
%  --- Output ---
% F: N x T matrix where N is the number of ROI, T is the number of frames in the Tiff files. 
% 
% example)
% 
% load('C:\home\ImagingExperiments\B6J647\F_B6J647_CREBRC107_20170902_d2_plane1_Ch1_Nk26.mat');
% Ipix = {stat.ipix};
% folder = get_tiff_location_from_ops(ops,'reg2P_kh004');
% folder(1)='D'; % change to D drive. Comment out if necessary.
% S=dir(fullfile(folder,'*.tif'));
% TiffNames= {S.names}; % assuming that Tiffs are correctly ordered by their numbers.
% 
% FullTiffNames=fullfile(folder,TiffNames);
% x_range = ops.xrange;
% y_range = ops.yrange;
% 
% F = get_F_mlttiff(Ipix, TiffNames,x_range,y_range);
%% by KH 20171010

Nroi = numel(Ipix); 
mt=mlttiff(FullTiffNames); 

if nargin >=4
    x_range = varargin{1};
    y_range = varargin{2};
else
    x_range = 1:mt.width;
    y_range = 1:mt.height;
end
%%
tic
F = nan(Nroi, mt.TotalNSeries , 'single');
for ii=1:length(mt.nSeries)
    index = mt.CumNSeries(ii)+[1:mt.nSeries(ii)];
    data = mt.get_stacks(index);
    data = data(y_range,x_range,:); 
    data = reshape(data,length(x_range)*length(y_range),[]);
    
    for kk = 1:Nroi
       ipix = Ipix{kk} ;
       if ~isempty(ipix)      
            F(kk,index) = mean(data(ipix,:), 1);       
       end
    end    
    fprintf('Frame %d/%d done in time %2.2f \n', index(end),mt.TotalNSeries, toc);
end



