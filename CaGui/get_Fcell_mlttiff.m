function F = get_Fcell_mlttiff(Ipix, FullTiffNames,x_range,y_range,Nframes)
%  A wrapper function to get Fcell from multiple tiff files.
% 
%  Note that, Fcell is different from F. 
%  In Fcell, there is a time gap between Fcell{1} and Fcell{2} during imaging.
%  F is concatenated version of Fcell, which is useful to calculate mean, std. 
% 
%  Usage1) Fcell = get_Fcell_mlttiff(Ipix, FullTiffNames,x_range,y_range,Nframes)
% 
%  --- Input1 ---
%  Ipix: 1xN or Nx1 cell array which contains indices of N- ROI. 
%  FullTiffNames: a cell of multiple tiff file names including all the path
%  x_range, y_range: Ipix is (often the case) indices in registered movie, 
%  which is smaller than the original Tiff file. To define the boundary of
%  the registered image, define x_range and y_range as 1 x W, 1xH vector.
%  If they are empty, x_range is assume to be 1:Width, and y_range = 1:Height of tiff file. 
% 
%  Nframes: a vector indicating the continuous frames without gap. 
%  for example, [30000, 30000] means there is a temporal gap between 30000
%  and 30001 frames, thus Fcell{1} contains 1 to 30000, and Fcell{2}
%  contains 30001 to 60000 frames.
% 
%  --- Output ---
% Fcell: 1 x length(Nframes) cell. Each cell contains N x Nframes(i) matrix
% 
% Example)
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
% Fcell = get_Fcell_mlttiff(Ipix, TiffNames,x_range,y_range,ops.Nframes);
%% by KH 20171010

F = get_signals_mlttiff(Ipix,FullTiffNames,x_range,y_range);

csumNframes = [0 cumsum(Nframes)];
Fcell = cell(1, length(Nframes));
for ii = 1:length(Nframes)
    Fcell{ii} = F(:, csumNframes(ii) + (1:Nframes(ii)));
end


