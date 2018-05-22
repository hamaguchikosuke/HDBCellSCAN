function [db,ops0,clustrules] = make_db_B6N792()
%% First, set the Root storage of your data
% ops0.RootStorage   = 'L:\home\ImagingData\Kosuke\DirectedLick\';
% ops0.RootStorage   = 'D:\home\ImagingData\Kosuke\DirectedLick';
ops0.RootStorage   = 'C:\home\GitHub\HDBCellScan\HDBCellScan\TestData';
ops0.useGPU                 = 0; % if you can use an Nvidia GPU in matlab this accelerates registration approx 3 times. You only need the Nvidia drivers installed (not CUDA).

%% clustrules.MinClust defines the size of cluster. Here I used pi*radius^2
clustrules.Compact          = 3; %  some cells with large dendrites could have values >2. I turned off filtering process in 
clustrules.diameter         = 12; % expected diameter of cells (used for 0.25 * pi/4*diam^2 < npixels < 30*pi/4*diam^2)
clustrules.npix_fraclow             = getOr(clustrules, {'npix_fraclow'}, (8/clustrules.diameter)^2);
clustrules.npix_frachigh            = getOr(clustrules, {'npix_frachigh'}, 50); %

% basic parameters for HDBSCAN based ROI detection 
clustrules.MinClust         = round(pi*clustrules.diameter^2/4); % to be used in HDBSCAN. 


%% ==== parameters in each database (db) ==== 
% Ca_channel: ROI will be detected in this channels. 
% Allign_channel: reference for movement correction. 
% Other channels will be registered referencing this channels movement data

%  Please put the data in the following structure. 
% RootStorage\mouse_name\date\expts1
% RootStorage\mouse_name\date\expts2
% 
% If a same plane is imaged but with any time gap, put Tiff files in different expts folder. 
% If a same plane is imaged but without any time gap, put Tiff files in the same expts folder.
% In db structure, movie in the sessions can be put together like
% db(i).expts = {'expts1','expts2'}
% Then registered images will be saved under 
% RootStorage\mouse_name\date\expts1_expts2

i = 0;


i = i+1;
db(i).mouse_name    = 'B6N792';
db(i).date          = '20180503';
% db(i).expts         = 'Test';
db(i).expts         = {'1_530um_2_530um'};
db(i).Ca_channel    = [1];
db(i).Align_channel   = 1; 
db(i).nchannels     = 1;
db(i).nplanes       = 1; 
db(i).comments      = 'Day85, second day of random session. ALM deep';
db(i).RegDirs       = {}; % blank it for the first time. if there are already registered files, set the parent directory here
db(i).doRegistration = 1; % set to 1 for the first time of analysis.
db(i).BFile = '';
db(i).RoughMeanTimeInBFile = NaN; % Give me one time point when the scanning was running during the behavior file. 


% 
% For 30FPS imaging
ops0.nFramesAvgForSVD       = 25; 
% For multiplane imaging, 5FPS
% ops0.nFramesAvgForSVD       = 6; 



% example extra entries
% db(i).AlignToRedChannel= 1;
% db(i).BiDiPhase        = 0; % adjust the relative phase of consecutive lines
% db(i).nSVD             = 1000; % will overwrite the default, only for this dataset
% db(i).comments      = 'this was an adaptation experiment';
% db(i).expred        = [4]; % say one block which had a red channel 
% db(i).nchannels_red = 2; % how many channels did the red block have in total (assumes red is last)