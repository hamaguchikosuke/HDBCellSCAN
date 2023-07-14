function varargout = HDBCellScan_Master_v03(varargin)
% HDBCELLSCAN_MASTER_V03 MATLAB code for HDBCellScan_Master_v03.fig
%      HDBCELLSCAN_MASTER_V03, by itself, creates a new HDBCELLSCAN_MASTER_V03 or raises the existing
%      singleton*.
%
%      H = HDBCELLSCAN_MASTER_V03 returns the handle to a new HDBCELLSCAN_MASTER_V03 or the handle to
%      the existing singleton*.
%
%      HDBCELLSCAN_MASTER_V03('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HDBCELLSCAN_MASTER_V03.M with the given input arguments.
%
%      HDBCELLSCAN_MASTER_V03('Property','Value',...) creates a new HDBCELLSCAN_MASTER_V03 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before HDBCellScan_Master_v03_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to HDBCellScan_Master_v03_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help HDBCellScan_Master_v03

% Last Modified by GUIDE v2.5 17-May-2018 15:26:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @HDBCellScan_Master_v03_OpeningFcn, ...
                   'gui_OutputFcn',  @HDBCellScan_Master_v03_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before HDBCellScan_Master_v03 is made visible.
function HDBCellScan_Master_v03_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to HDBCellScan_Master_v03 (see VARARGIN)

% Choose default command line output for HDBCellScan_Master_v03
handles.output = hObject;
% 
% toolbox_path = 'C:\home\matlab_svn\Suite2P';
% if exist(toolbox_path, 'dir')
% 	addpath(genpath(toolbox_path)) % add local path to the toolbox
% else
% 	error('toolbox_path does not exist, please change toolbox_path');
% end
% 



% Update handles structure
guidata(hObject, handles);

% UIWAIT makes HDBCellScan_Master_v03 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = HDBCellScan_Master_v03_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pb_loadDB.
function pb_loadDB_Callback(hObject, eventdata, handles)
% hObject    handle to pb_loadDB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try 
     [dbfile, dbpath]= uigetfile(fullfile(handles.DB.path,handles.DB.file),'Select a database m-file');
catch
    [dbfile, dbpath] = uigetfile({'*.m'},'Select a database m-file');
end
handles.DB.file = dbfile;
handles.DB.path = dbpath;

 addpath(dbpath);
 [~,dbfile,ext]=fileparts(dbfile);
 try 
     [db,ops,clustrules] = feval(dbfile);
 catch
     errordlg('Please update database file, it requires clustrules as third output');
 end
 ops.RootStorage=handles.DB.path; % overwrite the RootStorage,. 
 set(handles.text_dbname,'String',handles.DB.file);
 
 handles.DB.db = db;
 handles.DB.ops_original = ops;
 handles.DB.clustrules = clustrules;
 
 % check all the directory and returns the processed states.
 handles=pb_Scan_Callback(hObject, eventdata, handles);
%  set(handles.lb_sessions,'String',{db.date},'Value',1);

 guidata(hObject,handles);
 
% --- Executes on selection change in lb_sessions.
function lb_sessions_Callback(hObject, eventdata, handles)
% hObject    handle to lb_sessions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=popup_Methods_Callback(handles.popup_Methods, eventdata, handles);

% Hints: contents = cellstr(get(hObject,'String')) returns lb_sessions contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lb_sessions
contents = cellstr(get(hObject,'String')) ;
val = get(hObject,'Value');
% selected_name = contents{val};
% fprintf('===============================\n')
% fprintf('=== %s ===\n',contents{val});
% fprintf('===============================\n')
dbs=handles.DB.db(val);
for ii=1:length(dbs)
    fprintf('==================================================\n')
    fprintf('=== %s ===\n',contents{val(ii)});
   
    db=dbs(ii);    
    ops = handles.DB.ops_original;

    % scan directories to check the progress of data analysis.
    [ops,handles]=scan_db(ops,db,handles);
    
    handles.DB.ops = ops;
    handles=update_TxtStatus(hObject,eventdata,handles);
    fprintf('==================================================\n')
end
guidata(hObject,handles);


function [ops,handles]=scan_db(ops,db,handles)
% 1. Check whether image registration is done. 
% We can tell this if there are Tiff files under
% RootDir\AnimalID\SessionID\ContinuousImage\plane#_ch#\

ops = init_ops(ops);

ops=build_ops_for_HDBCellSCAN(db,ops);

% check registered tiff file exists
LatestDate=cell(length(ops.fsroot),2);

ops.process.RegTiffDone = true;
for ii=1:length(ops.fsroot)
    if isempty(ops.fsroot{ii})
        ops.process.RegTiffDone = false;  
        continue; 
    end
  % just check the first one
    tmp_fsroot = ops.fsroot{ii}(1);
    
    if ~exist(tmp_fsroot.name)
        ops.process.RegTiffDone = false;  
    else
        S=dir(tmp_fsroot.name);
%         if LatestDate{ii,1}<S.datenum, LatestDate(ii)
        if isfield(tmp_fsroot,'folder')
        ops.RegTiffPath{ii}=tmp_fsroot.folder; 
        else
            [folder,~,~]=fileparts(tmp_fsroot.name);
            ops.RegTiffPath{ii}=folder; 
        end
    end
end


[ops,handles]=update_PlaneChString(ops,handles);

 % 2. Check ROI is detected
% We can check this by inspecting whether 
% 1) there is a SVD file, SVDroi_<mouse_name>_<date>_plane<#>_ch<#>
%  which is used to detect ROIs, and then
% 2) there is ROI detection file,  ROI_<mouse_name>_<date>_plane<#>_ch<#>_Nk<ops.Nk0>.mat
% or  F_<mouse_name>_<date>_plane<#>_ch<#>_Nk<ops.Nk0>.mat (old suite2P)
% 

% update method  
handles=popup_Methods_Callback(handles.popup_Methods, [], handles);

switch handles.DB.ops.selected_analysis_ver 
    case 'ver20171103'
        ops.SVDSearchString = sprintf('SVDroi_%s_%s_plane*Ch*.mat',ops.mouse_name,ops.date);
        ROISearchStringBase = 'ROI_%s_%s_plane*ch*_MinClust%s.mat';
        ops.ROISearchString = sprintf(ROISearchStringBase,ops.mouse_name,ops.date,'*');
        ops.FSignalSearchString = sprintf('Fsig_%s_%s_plane*ch*_MinClust%s.mat',ops.mouse_name,ops.date,'*');
        ops.ProcSigSearchString = sprintf('Fsig_%s_%s_plane*ch*_MinClust%s_proc.mat',ops.mouse_name,ops.date,'*');
        ops.ProcSpkSearchString = sprintf('Fsig_%s_%s_plane*ch*_MinClust%s_procSpk.mat',ops.mouse_name,ops.date,'*');
    case 'Suite2P_2016'
        ops.SVDSearchString = sprintf('SVDroi_%s_%s_plane*.mat',ops.mouse_name,ops.date);
        ROISearchStringBase = 'F_%s_%s_plane*ch*_Nk%d.mat';
        ops.ROISearchString = sprintf(ROISearchStringBase,ops.mouse_name,ops.date,ops.Nk0);
        ops.FSignalSearchString = sprintf('F_%s_%s_plane*ch*_Nk%d.mat',ops.mouse_name,ops.date,ops.Nk);
        ops.ProcSigSearchString = sprintf('F_%s_%s_plane*ch*_Nk%d_proc.mat',ops.mouse_name,ops.date,ops.Nk);
        ops.ProcSpkSearchString = sprintf('F_%s_%s_plane*ch*_Nk%d_procSpk.mat',ops.mouse_name,ops.date,ops.Nk);
    otherwise
        error('Unknown Method %s',handles.DB.ops.selected_analysis_ver );
end

SVDFile = dir(fullfile(ops.ResultsSavePath,ops.SVDSearchString));
if ~isempty(SVDFile),   ops.SVDFile = SVDFile.name;  ops.process.SVDDone = true; ops.process.SVDdatenum=SVDFile.datenum;
else                    ops.SVDFile =  [];           ops.process.SVDDone = false;ops.process.SVDdatenum=0;
end

ROIFile = dir(fullfile(ops.ResultsSavePath,ops.ROISearchString));
if ~isempty(ROIFile),   ops.ROIFile = ROIFile.name;  ops.process.ROIDone = true; ops.process.ROIdatenum=ROIFile.datenum;
else                    ops.ROIFile =  [];           ops.process.ROIDone = false;ops.process.ROIdatenum=0;
end
  
% 3. Check signal is calculated
% We can check this by inspecting whether 
% there is F file,  Fsig_<mouse_name>_<date>_plane<#>_ch<#>_MinClust<#>.mat

FSignalFile = dir(fullfile(ops.ResultsSavePath,ops.FSignalSearchString));
if ~isempty(FSignalFile),   ops.FSignalFile = FSignalFile.name;  ops.process.FsignalDone = true;  ops.process.Fsignaldatenum = FSignalFile.datenum; 
else                        ops.FSignalFile =  [];               ops.process.FsignalDone = false; ops.process.Fsignaldatenum = 0; 
end
 

% 4. Check manual ROI selection is done. 
% We can check this by inspecting whether 
% there is F file,  Fsig_<mouse_name>_<date>_plane<#>_ch<#>_MinClust<#>_proc.mat

ProcFile = dir(fullfile(ops.ResultsSavePath,ops.ProcSigSearchString));
if ~isempty(ProcFile),   ops.ProcFile = ProcFile.name;  ops.process.ProcSignalDone = true; ops.process.ProcSignaldatenum=ProcFile.datenum;
else                     ops.ProcFile =  [];            ops.process.ProcSignalDone = false;ops.process.ProcSignaldatenum=0;
end
   
% 5. Check Fsig is merged with behavior file 
% We can check this by inspecting whether 
% there is F file,  Fsig_<mouse_name>_<date>_plane<#>_ch<#>_MinClust<#>_procSpk.mat

ProcSpkFile = dir(fullfile(ops.ResultsSavePath,ops.ProcSpkSearchString));
if ~isempty(ProcSpkFile),   ops.ProcSpkFile = ProcSpkFile.name;  ops.process.ProcSpkDone = true; ops.process.ProcSpkdatenum = ProcSpkFile.datenum;
else                        ops.ProcSpkFile =  [];               ops.process.ProcSpkDone = false;ops.process.ProcSpkdatenum = 0;
end
   
%% 
ops.ScanFileNames={'SVD','ROI','Fsignal','ProcSignal','ProcSpk'};
tmptb=struct2table(ops.process);
ScanFiledatenums=cellfun(@(x) sprintf('%sdatenum',x),ops.ScanFileNames,'UniformOutput',false);
[latestdatenum,latestInd]=max([tmptb(:,ScanFiledatenums).Variables]);
ops.process.latestInd    =latestInd;
ops.process.latestdatenum=latestdatenum;
ops.process.latestFile    =ops.ScanFileNames{latestInd};
% non-exist file is 0. exit file is 1. set the latest one as 2. This will be used to select the label in HDBCellSCAN_master list. 
ops.process.(sprintf('%sDone',ops.process.latestFile))=ops.process.(sprintf('%sDone',ops.process.latestFile))+1;

if ops.process.RegTiffDone, fprintf('1.Registered Image:\tFound in\t%s\n',ops.RegTiffPath{:}), else fprintf('1.Registered Image:\tNot Found\n'); end
if ops.process.SVDDone, fprintf('                  (SVD File:\tFound in\t%s)\n',ops.SVDFile), else fprintf('SVD File:\tNot Found\n'); end
if ops.process.ROIDone, fprintf('2.ROI File:\tFound in\t%s\n',ops.ROIFile), else fprintf('2.ROI File:\tNot Found\n'); end
if ops.process.FsignalDone, fprintf('3.Signal File:\tFound in\t%s\n',ops.FSignalFile), else fprintf('3.Fsignal File:\tNot Found\n'); end
if ops.process.ProcSignalDone, fprintf('4.Proc File:\tFound in\t%s\n',ops.ProcFile), else fprintf('4.Proc File:\tNot Found\n'); end
if ops.process.ProcSpkDone, fprintf('5.ProcSpk File:\tFound in\t%s\n',ops.ProcSpkFile), else fprintf('5.ProcSpk File:\tNot Found\n'); end


function [ops,handles]=update_PlaneChString(ops,handles)
% check fast reg file folder is created
handles=popup_Methods_Callback(handles.popup_Methods, [], handles);
switch handles.DB.ops.selected_analysis_ver 
    case 'ver20171103'
        S=dir(fullfile(ops.ResultsSavePath,'plane*ch*'));
    case 'Suite2P_2016'
        S=dir(fullfile(ops.ResultsSavePath,'plane*'));
     otherwise
        error('Unknown Method %s',handles.DB.ops.selected_analysis_ver );
end


if isempty(S)
    ops.process.RegTiffDone = false;
else
    DirInd = find([S.isdir]);
    PlaneCh=S(DirInd);
end

Val=get(handles.lb_PlaneCh,'Value');
    if isempty(Val)
        Val=1:length(S);
    end
    
% Before image registration, plane_ch dir does not exist. After collecting processed data into a single hard disk, 
% tiff files (except for x4 stacked images) are not copied, and this
% program will take it as "Not registered". In this case, existence of
% plane_ch folder will tell the user that image registered data alredy
% exists. 

if ops.process.RegTiffDone || ~isempty(S) 
    set(handles.lb_PlaneCh,'String',{PlaneCh.name},'Value',Val );
%     set(handles.lb_PlaneCh,'String',{PlaneCh.name},'Value', 1:length(S) );
  
end

function handles=update_TxtStatus(hObject,eventdata,handles)

ops = handles.DB.ops;
if ops.process.RegTiffDone,      set(handles.text_ImageReg,'String','Found');
else                             set(handles.text_ImageReg,'String','Not Found'); end


if ops.process.ROIDone, set(handles.text_GetROI,'String','Found'), 
else                    set(handles.text_GetROI,'String','Not Found'),  end

% if ops.process.ROIDone, fprintf('ROI File:\tFound in\t%s\n',ops.ROIFile), else fprintf('ROI File:\tNot Found\n'); end

if ops.process.FsignalDone, set(handles.text_GetSignal,'String','Found'), 
else,                       set(handles.text_GetSignal,'String','Not Found'); end

if ops.process.ProcSignalDone, set(handles.text_ProcFile,'String','Found'), 
else,                       set(handles.text_ProcFile,'String','Not Found'); end

if ops.process.ProcSpkDone, set(handles.text_ProcSpk,'String','Found'), 
else,                       set(handles.text_ProcSpk,'String','Not Found'); end

% We may also need to update PlaneCh panel for each data entry. 
% Check whether this does not cause any problem .
[ops,handles]=update_PlaneChString(ops,handles);

% --- Executes during object creation, after setting all properties.
function lb_sessions_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lb_sessions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in rb_imagereg.
function rb_imagereg_Callback(hObject, eventdata, handles)
% hObject    handle to rb_imagereg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rb_imagereg


% --- Executes on button press in rb_GetRoi.
function rb_GetRoi_Callback(hObject, eventdata, handles)
% hObject    handle to rb_GetRoi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rb_GetRoi


% --- Executes on button press in rb_GetSignal.
function rb_GetSignal_Callback(hObject, eventdata, handles)
% hObject    handle to rb_GetSignal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rb_GetSignal



% --- Executes on button press in pb_Run.
function pb_Run_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if numel(handles.lb_sessions,'Value')>1
ButtonName = questdlg({'==== This is a batch analysis mode ====',...
    'You have two options.',...
    'Overwrite: redo all the analysis and overwrite.', ...
    'Keep : skip if the analysis is done already'}, ...
    'Question', ...
    'Overwrite', 'Keep', 'Cancel', 'Keep');
else
    % in a single run mode, always follow the users instruction.
    handles.DB.redo_analysis=1;
end

switch ButtonName,
    case 'Overwrite'
        disp('Overwrite: Re-do all the analysis.');
        handles.DB.redo_analysis=1;
    case 'Keep',
        disp('Keep: Analyze unprocessed data.');
        handles.DB.redo_analysis=0;
    case 'Cancel',
        disp('User canceled.')
        return;
end % switch

contents = get(handles.lb_sessions,'String');   
val = get(handles.lb_sessions,'Value');   
dbs = handles.DB.db(val);

for ii=1:length(dbs)
    fprintf('==================================================\n')
    fprintf('=== %s ===\n',contents{val(ii)});
   
    db=dbs(ii);    
    ops = handles.DB.ops_original;

    % scan directories to check the progress of data analysis.
    [ops,handles]=scan_db(ops,db,handles);
    handles.DB.ops = ops;
    SingleRun(handles.pb_Run,eventdata,handles,db);
  
    handles=update_TxtStatus(hObject,eventdata,handles);
    fprintf('==================================================\n')
end

fprintf('============================== Complete ==============================\n');

function SingleRun(hObject,eventdata,handles,db)
TSTART = tic;

% db = handles.DB.db;
ops0 = update_ops(handles.DB.ops_original,handles.DB.ops);  
clustrules = handles.DB.clustrules;

ops0.doRegistration  = get(handles.rb_imagereg,'Value');
ops0.getROIs         = get(handles.rb_GetRoi,  'Value');
ops0.RegFile_xtime   = getOr(db, {'RegFile_xtime'},4); % number of tiffs to average when writing time averaged registered tiffs (slow)
ops0.DoGetSignal     = get(handles.rb_GetSignal,'Value');
ops0.DoUpdateSignal  = get(handles.rb_FinalizeSignal,'Value');
%%%%%%%%%%%% image registration and svd computation %%%%%%%%%%%%
update_reg(ops0,handles);
% registered tiff are stored in 
% RegFileTiffLocation\mouse_name\session(date)\expts\plane#_ch#\
% SVD file is in SVDroi_<mouse_name>_<date>_plane<#>_ch<#>.mat
  

% After registration, Plane#_Ch# folder is created. Update it. 
% [~,handles]=update_PlaneChString(ops0,handles);
lb_sessions_Callback(handles.lb_sessions, eventdata, handles);

%%%%%%%%% get ROIs for selected plane, view, and channel %%%%%%%%
% by using SVD file computed above. 
update_roi(ops0,handles)
% ROI file is stored in  
% ROI_<mouse_name>_<date>_plane<#>_ch<#>.mat
     

%%%%%%%%% get signals for selected plane, view, channel %%%%%%%%
% by using ROI file and registered Tiff computed above, 
calc_signal(ops0,handles);
% Signal file is stored in
% Fsig_<mouse_name>_<date>_plane<#>_ch<#>.mat

%%%%%%%%% manual check  ROI %%%%%%%%
% call Roi curating program (RoiGui)
manual_roi(ops0,handles);


%%%%%%%%% get signals for selected plane, view, channel %%%%%%%%
% by using Proc file and registered Tiff computed above, 
update_signal(ops0,handles);
% Proc file is stored in
% Fsig_<mouse_name>_<date>_plane<#>_ch<#>_proc.mat

toc(TSTART);

%%%%% cleanup %%%%%
try
    for ii=1:length(ops1)       
        if ops1{ii}.DeleteBin
            fclose('all');
            delete(ops1{ii}.RegFile);        % delete temporary bin file
        end
    end
end
gong(5,400,1);pause(1);gong(5,400,1);pause(1);gong(5,400,1);

function ops1=update_reg(ops0,handles)
 
%%%%%%%%%%%% image registration and svd computation %%%%%%%%%%%%
ops1 = [];
if get(handles.rb_imagereg,'Value') && ...
        (handles.DB.redo_analysis || (~ops0.process.RegTiffDone || ~ops0.process.SVDDone) )
        ops1         = reg2P_kh007(ops0);    %  ver.7. GrinLens mode is added. 
%      ops1         = reg2P_kh006(ops0);  %  ver.6. GrinLens mode is added. 
%       ops1         = reg2P_kh005(ops0);  % ver.5
     % At this point, ops1 contains Plane x View x Channel size. 
    % get_SVD function relies on temporally generated temp_file, so always use right after reg2P. 
    svd_start=tic; 
    ops1         = get_svd_20171102(ops1); % get svd 
    fprintf('SVD computation time %f\n',toc(svd_start));
end

function ops1 = update_roi(ops0,handles)
%%%%%%%%% get ROIs for selected plane, view, and channel %%%%%%%%

% update the following: handles.DB.ops.selected_analysis_ver 
handles=popup_Methods_Callback(handles.popup_Methods, [], handles);

    
PlaneChString = get(handles.lb_PlaneCh,'String');
PlaneChString =  PlaneChString(get(handles.lb_PlaneCh,'Value'));

           
if  get(handles.rb_GetRoi,'Value') && ...
        (handles.DB.redo_analysis || ~ops0.process.ROIDone)
    
    ops1 = [];

    for pp=1:length(PlaneChString)
        SVDFile = sprintf('%s/SVDroi_%s_%s_%s.mat', ops0.ResultsSavePath, ...
            ops0.mouse_name, ops0.date, PlaneChString{pp});
        fprintf('Loading svd file: %s\n', SVDFile);
      
        if exist(SVDFile,'file'),        tmp=   load(SVDFile,  'ops','U','Sv');
        else
            error('SVD file %s NOT FOUND.',SVDFile);
        end
        
       % tmp.ops contains previous options. ops0 inherits latest options.
       % ops1{pp} needs to be constructed based on ops0, but otherwise use old values 
%        ops0.PlaneChString= PlaneChString{pp};
        % just to make sure PlaneID is updated. 
       try ops0 = rmfield(ops0,{'PlaneID','ChannelID'});end
       ops1{pp} = update_ops(ops0,tmp.ops);    
       % The following two are exceptions, they will overwrite the
       % parameter 
       ops1{pp}.ChannelID=tmp.ops.ChannelID;
       ops1{pp}.PlaneID=tmp.ops.PlaneID;
       roi_start =tic;
       switch handles.DB.ops.selected_analysis_ver
           case 'ver20171103'
              ops1{pp} = get_roi_20171102(ops1{pp},tmp.U,tmp.Sv,handles.DB.clustrules);
           case 'Suite2P_2016'
               ops_tmp=ops1{pp};
               ops_tmp.diameter = getOr(ops_tmp,'diameter',handles.DB.clustrules.diameter);
               ops_tmp.iplane = getOr(ops_tmp,'iplane',1);
               ops_tmp.Nk0=500;
               ops_tmp.Nk=250;
               ops_tmp.clustModel='neuropil';
               ops1{pp} = get_roi_20171102(ops_tmp,tmp.U,tmp.Sv,handles.DB.clustrules);
%               [ops1{pp}, stat, res]  = fast_clustering_with_neuropil(ops_tmp,tmp.U,tmp.S);
              
           otherwise
               error('Unknown Method %s',handles.DB.ops.selected_analysis_ver );
       end
      
       fprintf('ROI computation time: %fsec\n',toc(roi_start));
    end

end

function ops1 = calc_signal(ops0,handles)
%%%%%%%%% get signals for selected plane, view, channel %%%%%%%%
% by using ROI file and registered Tiff computed above two. 
% update the following: handles.DB.ops.selected_analysis_ver 
handles=popup_Methods_Callback(handles.popup_Methods, [], handles);

PlaneChString = get(handles.lb_PlaneCh,'String');
PlaneChString =  PlaneChString(get(handles.lb_PlaneCh,'Value'));

    
if  get(handles.rb_GetSignal,'Value') && ...
        (handles.DB.redo_analysis || ~ops0.process.FsignalDone)
        
       ops1 = [];
  for pp=1:length(PlaneChString)
      
      switch handles.DB.ops.selected_analysis_ver
          case 'ver20171103'
              MinClustSize = handles.DB.clustrules.MinClust;
              ROIFile = sprintf('%s/ROI_%s_%s_%s_MinClust%d.mat', ops0.ResultsSavePath, ...
                  ops0.mouse_name, ops0.date, PlaneChString{pp},MinClustSize);
          case 'Suite2P_2016'
                ops0.iplane = getOr(ops0,'iplane',1);
              ROIFile=sprintf('%s/F_%s_%s_plane%d_Ch1_Nk%d.mat', ops0.ResultsSavePath, ...
                  ops0.mouse_name, ops0.date, ops0.iplane, ops0.Nk);
          otherwise
              error('Unknown option');
      end
        
        if exist(ROIFile,'file')
            fprintf('Loading ROI file: %s\n', ROIFile);
            tmp=load(ROIFile,'ops', 'res', 'stat', 'stat0', 'res0', 'clustrules');
            
            
            % tmp.ops contains previous options. ops0 inherits latest options.
            % ops1{pp} needs to be constructed based on ops0, but otherwise use old values
            
            % just to make sure PlaneID is updated. 
%             ops0 = rmfield(ops0,'PlaneID');
            ops1{pp} = update_ops(ops0,tmp.ops);
            
            ops1{pp}.PlaneChString = PlaneChString{pp};
            
            neuropilSub    = 'surround'; %getOr(ops1, {'neuropilSub'}, 'surround')
            
            switch neuropilSub
                case 'surround'
                    get_signals_and_neuropil_kh(ops1{pp}, ops1{pp}.PlaneChString);
                otherwise
                    error('Unknown option!');
            end
            
            
        else
            errordlg(sprintf('%s does not exist.', ROIFile));
        end
        
                   
  end 
  
end
fprintf('==== Done ====\n');

function manual_roi(ops0,handles)

%%%%%%%%% get ROIs for selected plane, view, and channel %%%%%%%%
PlaneChString = get(handles.lb_PlaneCh,'String');
PlaneChString =  PlaneChString(get(handles.lb_PlaneCh,'Value'));

           
if  get(handles.rb_ManualROI,'Value') && ...
        (handles.DB.redo_analysis || ~ops0.process.ProcSignalDone)
    
    ops1 = [];

    for pp=1:length(PlaneChString)
        MinClustSize = handles.DB.clustrules.MinClust;
        FsigFile = sprintf('%s/Fsig_%s_%s_%s_MinClust%d.mat', ops0.ResultsSavePath, ...
            ops0.mouse_name, ops0.date, PlaneChString{pp},MinClustSize);
        
        
        if exist(FsigFile,'file')
            fprintf('Please load Fsig file: %s\n', FsigFile);
%             pb_LoadPlane_Callback(hObject, eventdata, h,varargin)
            h=guidata(RoiGui_007);
%             RoiGui_007('pb_LoadPlane_Callback',h.figure1,[],h,FsigFile);
            RoiGui_008('pb_LoadPlane_Callback',h.figure1,[],h,FsigFile);
            % <To Dos>: load Fsig file automatically.
            
        else
            errordlg(sprintf('%s does not exist.', FsigFile));
        end
    end

end


function ops1 = update_signal(ops0,handles)
%%%%%%%%% get signals for selected plane, view, channel %%%%%%%%
% by using Proc file and registered Tiff computed above two. 

PlaneChString = get(handles.lb_PlaneCh,'String');
PlaneChString =  PlaneChString(get(handles.lb_PlaneCh,'Value'));

    
if  get(handles.rb_FinalizeSignal,'Value') && ...
        (handles.DB.redo_analysis || ~ops0.process.ProcSpkDone)
    ops1 = []; 
           
       
  for pp=1:length(PlaneChString)
      MinClustSize = handles.DB.clustrules.MinClust;
      
      ProcFile = fullfile(ops0.ResultsSavePath,...
          sprintf('Fsig_%s_%s_%s_MinClust%d_proc.mat',ops0.mouse_name, ops0.date, PlaneChString{pp},MinClustSize));
      ProcSpkFile = fullfile(ops0.ResultsSavePath,...
          sprintf('Fsig_%s_%s_%s_MinClust%d_procSpk.mat', ops0.mouse_name, ops0.date, PlaneChString{pp},MinClustSize));
     
     % This ProcFile could depend on older MinClust parameter, in that case, failed to find the current ProcFile.
     % Use ops0.ProcFile 
      
     if ~exist(ProcFile,'file') && exist(fullfile(ops0.ResultsSavePath,ops0.ProcFile),'file')
         ProcFile=fullfile(ops0.ResultsSavePath,ops0.ProcFile);
     end
         
          
      if exist(ProcFile,'file') & handles.rb_updateSpikeOnly.Value
       
          % In case update (ProcSpkFile already exists but overwrite)
          proc2procspk_02(ProcFile);
      elseif exist(ProcFile,'file')
          fprintf('Loading Proc file: %s\n', ProcFile);
          tmp=load(ProcFile,'ops');

          % tmp.ops contains previous options. ops0 inherits latest options.
          % ops1{pp} needs to be constructed based on ops0, but otherwise use old values
           ops1{pp} = update_ops(ops0,tmp.ops);
          
% the following code ask user to decide where to save. When commented out, always saved under ops1{pp}.ResultsSavePath            
%           if ~strncmp(ops0.ResultsSavePath , ops1.ResultsSavePath,length(ops0.ResultsSavePath))
%               ButtonName = questdlg('ResultsSavePath conflict: likely you worked on a local _proc file, and update the signal by using the original data source. Which one to use?', ...
%                   'Conflict ', ...
%                   ProcFile, ops0.ResultsSavePath,'Cancel');
%               switch ButtonName
%                   case {ProcFile,ops0.ResultsSavePath}
%                       fprintf('Use %s\n',ButtonName);
%                   case 'Cancel'
%                       disp('User cancelled');
%                       return;
%               end % switch
%               
%           end
         
          ops1{pp}.PlaneChString = PlaneChString{pp};
          
          neuropilSub    = 'surround'; %getOr(ops1, {'neuropilSub'}, 'surround')
          
          switch neuropilSub
              case 'surround'
                  get_signals_and_neuropil_kh(ops1{pp}, ops1{pp}.PlaneChString);
              otherwise
                  error('Unknown option!');
          end
          
          % then, add deconvoluted data.
          proc2procspk_02(ProcFile);
          
          
      else
          errordlg(sprintf('%s does not exist.', ProcFile));
      end
      
  end 
end
fprintf('==== Done ====\n');

% --- Executes on button press in pb_Cancel.
function pb_Cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in lb_PlaneCh.
function lb_PlaneCh_Callback(hObject, eventdata, handles)
% hObject    handle to lb_PlaneCh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lb_PlaneCh contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lb_PlaneCh


% --- Executes during object creation, after setting all properties.
function lb_PlaneCh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lb_PlaneCh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_Methods.
function handles=popup_Methods_Callback(hObject, eventdata, handles)
% hObject    handle to popup_Methods (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_Methods contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_Methods
contents = cellstr(get(hObject,'String')) ;
handles.DB.ops.selected_analysis_ver = contents{get(hObject,'Value')};

% --- Executes during object creation, after setting all properties.
function popup_Methods_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_Methods (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in rb_ManualROI.
function rb_ManualROI_Callback(hObject, eventdata, handles)
% hObject    handle to rb_ManualROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rb_ManualROI


% --- Executes on button press in rb_FinalizeSignal.
function rb_FinalizeSignal_Callback(hObject, eventdata, handles)
% hObject    handle to rb_FinalizeSignal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rb_FinalizeSignal


% --- Executes on button press in pb_Scan.
function handles=pb_Scan_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Scan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%update which version to use
handles=popup_Methods_Callback(handles.popup_Methods, eventdata, handles);

dbs=handles.DB.db;

for ii=1:length(dbs)
    db=dbs(ii);
    ops = handles.DB.ops_original;
    
    % scan directories to check the progress of data analysis.
    [ops,handles]=scan_db(ops,db,handles);
    ProcessOrNot  = {'-','*','o'}; 
    
    ProcessOrNotString = sprintf('%s%s%s%s%s',...
        ProcessOrNot{ops.process.RegTiffDone+1},...
        ProcessOrNot{ops.process.ROIDone+1},...
        ProcessOrNot{ops.process.FsignalDone+1},...
        ProcessOrNot{ops.process.ProcSignalDone+1},...
        ProcessOrNot{ops.process.ProcSpkDone+1}  );
    
    Session_String = [sprintf('%-20s:',handles.DB.db(ii).date),sprintf('\t%s',ProcessOrNotString)];
    handles.DB.db(ii).session_string = Session_String;
    handles.DB.db(ii).scan_results = ops;
end

dbs=handles.DB.db;
set(handles.lb_sessions,'String',{dbs.session_string},'Value',1,'FontName','FixedWidth');

 guidata(hObject,handles);


% --- Executes on button press in rb_BatchSession.
function rb_BatchSession_Callback(hObject, eventdata, handles)
% hObject    handle to rb_BatchSession (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rb_BatchSession


% --- Executes on button press in rb_updateSpikeOnly.
function rb_updateSpikeOnly_Callback(hObject, eventdata, handles)
% hObject    handle to rb_updateSpikeOnly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rb_updateSpikeOnly
