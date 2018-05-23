%% setup path to matlab and python 
addpath(genpath('C:\home\GitHub\HDBCellScan')) % path to HDBCellScan folder
addpath(genpath('C:\home\GitHub\suite2P')) % path to suite2P 

P=pyversion;
if isempty(P)
    errordlg('Python is not found. Install Python or set PATH so that Windows can find it');
end
% Before run this code, please check whether you have installed HDBSCAN
% toolbox through Anaconda. Then, please find folder that includes hdbscan,
% such as 

% python_path = 'C:\Users\hamag\.conda\envs\khtest\Lib\site-packages'; 
python_path = 'C:\Users\hammer\AppData\Local\conda\conda\envs\CaImaging\Lib\site-packages'
P = py.sys.path;
append(P,python_path);

%% First, setup database m-file 
edit('C:\home\GitHub\HDBCellScan\TestData\make_db_B6N792');
%% Load the database m-file from master control 
HDBCellScan_Master_v03; 

% ==== Brief Instructions====
% 
% ---- analysis ----
% (1): ImageReg 
%   image registration and compute svd matrix from the image
%   This part utilized suite2P codes. 
% (2): GetRoi 
%   This is the heart of HDBCellScan. 
%   It internally calls python toolbox HDBSCAN. 
% (3): GetSignal
%   initial step of signal detection (fluorescence).
% (4): Gui-based ROI selection 
%   Inspect the quality of automatically found ROIs. 
%   Teach the program which ROIs are cell bodies, and which are noises,
%   dendrites.
%   Split an ROI into two if necessary.
% (5): Finalize the signal (F) and estimate spikes using equivalent algorithm of fast-oopsi. 
% 
% First, press [Load DB] button to load make_db_B6N792.m located under TestData. 
% 
% Second, select the database entries you want to analyze. 
% In this test database, there is only one entry "20180503"
% 
% Third, select the first three (ImageReg, GetRoi, GetSignal) analysis and press [RUN].
% You can choose 
% Overwrite: re-do analysis fromthe beginning.
% Keep:      skip the analysis if previously analyzed file exist.
% Cancel:    Abort analysis.
% 
% Then, select fourth button (GUI ROI Selection), then it will open a GUI
% to semi-automatically select ROIs.  (Please refer to another manual for
% this GUI). At the end, save proc file. 
% 
% Finally, finalize the sinal by turn on fifth button [Finalize Signal] and Run. 
% Done! 
%  
% The final output is
% Fsig_<AnimalName>_<Date>_<plane#_ch#>_MinClust#_procSpk.mat

%% Data structure in procSpk
% 
% ===== F ====
% F contains all the fluorescent and extracted action potential data.
% F.FcellNeu: Neuropil subtracted Fluoresecent data.
% F.FcellNeu{1}(cl.selected,:) are selected ROIs fluorescence signal
% through HDBCellSCAN_GUI.
% F.Spk     : estimated action potential data
% 
% 
% 
% 
% 
