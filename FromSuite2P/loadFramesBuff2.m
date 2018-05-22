function [frames, headers] = loadFramesBuff2(tiff, firstIdx, lastIdx, stride, ...
    ind_i,ind_j,temp_file)
%loadFramesBuff2 Loads the frames of a Tiff file into an array (Y,X,T)
%   MOVIE = loadFrames(TIFF, [FIRST], [LAST], [STRIDE], []) loads
%   frames from the Tiff file specified by TIFF, which should be a filename
%   or an already open Tiff object. Optionallly FIRST, LAST and STRIDE
%   specify the range of frame indices to load.
% 
%  Option2: load a subregion of the image
%  MOVIE = loadFrames(TIFF, [FIRST], [LAST], [STRIDE], ind_i,ind_j)
%  returns Movie(ind_i,ind_j,[FIRST:STRIDE:LAST]);
%  by KH 20161222

TiffFileName = tiff;
[~,TiffFileName,~]=fileparts(TiffFileName);

% initChars = overfprintf(0, 'Loading TIFF frame ');
warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');

warningsBackOn = onCleanup(...
  @() warning('on', 'MATLAB:imagesci:tiffmexutils:libtiffWarning'));


if nargin < 2 || isempty(firstIdx)
  firstIdx = 1;
end

if nargin < 3 || isempty(lastIdx)
  lastIdx = nFrames(tiff);
end

if nargin < 4 || isempty(stride)
  stride = 1;
end

if nargin>=7
   copyfile(tiff,temp_file) 
   tiff = temp_file;
end

if ischar(tiff)
  tiff = Tiff(tiff, 'r');
  closeTiff = onCleanup(@() close(tiff));
end

if nargout > 1
  loadHeaders = true;
else
  loadHeaders = false;
end

w = tiff.getTag('ImageWidth');
h = tiff.getTag('ImageLength');

if nargin < 5 || isempty(ind_i)
  ind_i = 1:h;
end

if nargin < 6 || isempty(ind_j)
  ind_j = 1:w;
end

dataClass = class(read(tiff));
TotalFrameNums = ceil((lastIdx - firstIdx + 1)/stride);
frames = zeros(length(ind_i), length(ind_j), TotalFrameNums, dataClass);
if loadHeaders
  headers = cell(1, TotalFrameNums);
end

nMsgChars = 0;
setDirectory(tiff, firstIdx);
waitCnt = 0;
waitH=waitbar(waitCnt,'Loading Tiffs ( / )...');
for t = 1:TotalFrameNums
  if t/TotalFrameNums>=waitCnt+0.05
      waitCnt=waitCnt+0.1;
      waitbar(waitCnt,waitH,{strrep(TiffFileName,'_','\_'),sprintf('Loading Tiffs (%d/%d)',t,TotalFrameNums)});
    %nMsgChars = overfprintf(nMsgChars, '%i/%i', t, nFrames);
  end
  
  tmp=read(tiff);
  frames(:,:,t) = tmp(ind_i,ind_j);
  
  if loadHeaders
    headers{t} = getTag(tiff, 'ImageDescription');
  end
  
  if t < TotalFrameNums
    for i = 1:stride
      nextDirectory(tiff);
    end
  end
end

close(waitH);
%overfprintf(initChars + nMsgChars, '');

end

