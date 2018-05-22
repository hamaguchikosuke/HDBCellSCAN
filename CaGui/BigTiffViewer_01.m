%  large Tiff file sequence viewer;
filenames = {'20161227_1_b6j498_2P_plane1_x5_1.tif',...
                      '20161227_1_b6j498_2P_plane1_x5_2.tif',...
                      '20161227_1_b6j498_2P_plane1_x5_3.tif'};
                  
% filepath  =  fileparts(handles.dat.filename); 
filepath  = 'F:\home\ImagingData\DualLickMice\b6j498\20161227\1\Plane1';
full_filenames = fullfile(filepath,filenames);
NumSeries=cellfun(@nFrames, full_filenames);
CumNumSeries = [0 cumsum(NumSeries)]; % add zero for later purpose.
Loaded = 0;
%%
myfigure('Big Tiff Viewer'); clf;
% index = 5650;

for index=1:100:4000
    
FileInd = find(CumNumSeries>index,1,'first')-1;
LocalInd = index-CumNumSeries(FileInd);

if Loaded ~= FileInd
    
    firstIdx=1;
    lastIdx=NumSeries(FileInd);
    stride = 1;
    ind_i = [];
    ind_j = [];
    
    IMG=loadFramesBuff2(full_filenames{FileInd}, ...
        firstIdx , lastIdx, stride, ...
        ind_i,ind_j);
    Loaded = FileInd;
end

image(IMG(:,:,LocalInd));
title(num2str(index));
drawnow;
end