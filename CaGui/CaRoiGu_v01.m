function varargout = CaRoiGu_v01(varargin)
% CAROIGU_V01 MATLAB code for CaRoiGu_v01.fig
%      CAROIGU_V01, by itself, creates a new CAROIGU_V01 or raises the existing
%      singleton*.
%
%      H = CAROIGU_V01 returns the handle to a new CAROIGU_V01 or the handle to
%      the existing singleton*.
%
%      CAROIGU_V01('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CAROIGU_V01.M with the given input arguments.
%
%      CAROIGU_V01('Property','Value',...) creates a new CAROIGU_V01 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CaRoiGu_v01_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CaRoiGu_v01_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CaRoiGu_v01

% Last Modified by GUIDE v2.5 16-Nov-2017 15:25:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CaRoiGu_v01_OpeningFcn, ...
                   'gui_OutputFcn',  @CaRoiGu_v01_OutputFcn, ...
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


% --- Executes just before CaRoiGu_v01 is made visible.
function CaRoiGu_v01_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CaRoiGu_v01 (see VARARGIN)

% Choose default command line output for CaRoiGu_v01
handles.output = hObject;

toolbox_path = 'C:\home\matlab_svn\Suite2P';
if exist(toolbox_path, 'dir')
	addpath(genpath(toolbox_path)) % add local path to the toolbox
else
	error('toolbox_path does not exist, please change toolbox_path');
end

handles.toolbox_path = 'C:\home\matlab_svn\suite2P\';
if exist(handles.toolbox_path, 'dir')
	addpath(genpath(handles.toolbox_path)) % add local path to the toolbox
else
	error('toolbox_path does not exist, please change toolbox_path');
end



% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CaRoiGu_v01 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = CaRoiGu_v01_OutputFcn(hObject, eventdata, handles) 
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
 
 handles.DB.db = db;
 handles.DB.ops_original = ops;
 handles.DB.clustrules = clustrules;
 
 set(handles.lb_sessions,'String',{db.date},'Value',1);

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
selected_name = contents{val};
fprintf('===============================\n')
fprintf('=== %s selected===\n',selected_name);
fprintf('===============================\n')
db=handles.DB.db(val);
ops = handles.DB.ops_original;

% 1. Check whether image registration is done. 
% We can tell this if there are Tiff files under
% RootDir\AnimalID\SessionID\ContinuousImage\plane#_ch#\

ops = init_ops(ops);
ops=build_ops3(db,ops);

% check registered tiff file exists
ops.process.RegTiffDone = true;
for ii=1:length(ops.fsroot)
    if ~exist(ops.fsroot{ii}.name), 
        ops.process.RegTiffDone = false;  
    else
        if isfield(ops.fsroot{1},'folder')
        ops.RegTiffPath{ii}=ops.fsroot{ii}.folder; 
        else
            [folder,~,~]=fileparts(ops.fsroot{ii}.name);
            ops.RegTiffPath{ii}=folder; 
        end
    end
end

% check fast reg file folder is created
S=dir(fullfile(ops.ResultsSavePath,'plane*ch*'));
if isempty(S)
    ops.process.RegTiffDone = false;
else
    DirInd = find([S.isdir]);
    PlaneCh=S(DirInd);
end

% ops.SearchTiffString=get_tiff_location_from_ops_20171102(ops);
% ops.FastRegTiffPath = fullfile(ops.ResultsSavePath, {PlaneCh.name});

% if isempty(ops.RegTiffPath),    ops.process.RegTiffDone = false;  % yet 
% else                            ops.process.RegTiffDone = true;   end


if ops.process.RegTiffDone
    set(handles.lb_PlaneCh,'String',{PlaneCh.name},'Value',1  );
end
% 2. Check ROI is detected
% We can check this by inspecting whether 
% 1) there is a SVD file, SVDroi_<mouse_name>_<date>_plane<#>_ch<#>
%  which is used to detect ROIs, and then
% 2) there is ROI detection file,  ROI_<mouse_name>_<date>_plane<#>_ch<#>_Nk<ops.Nk0>.mat
% or  F_<mouse_name>_<date>_plane<#>_ch<#>_Nk<ops.Nk0>.mat (old suite2P)
% 
% 

ops.SVDSearchString = sprintf('SVDroi_%s_%s_plane*Ch*.mat',ops.mouse_name,ops.date);
SVDFile = dir(fullfile(ops.ResultsSavePath,ops.SVDSearchString));
if ~isempty(SVDFile),   ops.SVDFile = SVDFile.name;  ops.process.SVDDone = true;
else                    ops.SVDFile =  [];           ops.process.SVDDone = false; end

switch handles.DB.ops.selected_analysis_ver 
    case 'ver20171103'
        ROISearchStringBase = 'ROI_%s_%s_plane*ch*_MinClust%s.mat';
        ops.ROISearchString = sprintf(ROISearchStringBase,ops.mouse_name,ops.date,'*');
    case 'Suite2P_2016'
        ROISearchStringBase = 'F_%s_%s_plane*ch*_Nk%d.mat';
        ops.ROISearchString = sprintf(ROISearchStringBase,ops.mouse_name,ops.date,ops.Nk0);
    otherwise
        error('Unknown Method %s',handles.DB.ops.selected_analysis_ver );
end

ROIFile = dir(fullfile(ops.ResultsSavePath,ops.ROISearchString));
if ~isempty(ROIFile),   ops.ROIFile = ROIFile.name;  ops.process.ROIDone = true;
else                    ops.ROIFile =  [];           ops.process.ROIDone = false; end
  
% 3. Check signal is calculated
% We can check this by inspecting whether 
% there is F file,  F_<mouse_name>_<date>_plane<#>_ch<#>_Nk<ops.Nk>.mat

ops.FSignalSearchString = sprintf('F_%s_%s_plane*ch*_Nk%d.mat',ops.mouse_name,ops.date,ops.Nk);
FSignalFile = dir(fullfile(ops.ResultsSavePath,ops.FSignalSearchString));
if ~isempty(FSignalFile),   ops.FSignalFile = FSignalFile.name;  ops.process.FsignalDone = true;
else                        ops.FSignalFile =  [];               ops.process.FsignalDone = false; end
 
if ops.process.RegTiffDone, fprintf('Registered Image:\tFound in\t%s\n',ops.RegTiffPath{:}), else fprintf('Registered Image:\tNot Found\n'); end
if ops.process.SVDDone, fprintf('SVD File:\tFound in\t%s\n',ops.SVDFile), else fprintf('SVD File:\tNot Found\n'); end
if ops.process.ROIDone, fprintf('ROI File:\tFound in\t%s\n',ops.ROIFile), else fprintf('ROI File:\tNot Found\n'); end
if ops.process.FsignalDone, fprintf('Signal File:\tFound in\t%s\n',ops.FSignalFile), else fprintf('Fsignal File:\tNot Found\n'); end

handles.DB.ops = ops;
handles=update_TxtStatus(hObject,eventdata,handles);
guidata(hObject,handles);

function handles=update_TxtStatus(hObject,eventdata,handles)

ops = handles.DB.ops;
if ops.process.RegTiffDone,      set(handles.text_ImageReg,'String','Found');
else                             set(handles.text_ImageReg,'String','Not Found'); end


if ops.process.ROIDone, set(handles.text_GetROI,'String','Found'), 
else                    set(handles.text_GetROI,'String','Not Found'),  end

% if ops.process.ROIDone, fprintf('ROI File:\tFound in\t%s\n',ops.ROIFile), else fprintf('ROI File:\tNot Found\n'); end

if ops.process.FsignalDone, set(handles.text_GetSignal,'String','Found'), 
else,                       set(handles.text_GetSignal,'String','Not Found'); end



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

DoReg =         get(handles.rb_imagereg,'Value');
DoCalcROI =     get(handles.rb_GetRoi,  'Value');
DoGetSignal =   get(handles.rb_GetSignal,'Value');

db = handles.DB.db;
ops0 = update_ops(handles.DB.ops_original,handles.DB.ops);  
clustrules = handles.DB.clustrules;

ops0.doRegistration  = DoReg;
ops0.getROIs         = DoCalcROI;
ops0.RegFile_xtime   = getOr(db, {'RegFile_xtime'},4); % number of tiffs to average when writing time averaged registered tiffs (slow)
ops0.DoGetSignal     = DoGetSignal;

%%%%%%%%%%%% image registration and svd computation %%%%%%%%%%%%
update_reg(ops0);
% registered tiff are stored in 
% RegFileTiffLocation\mouse_name\session(date)\expts\plane#_ch#\
% SVD file is in SVDroi_<mouse_name>_<date>_plane<#>_ch<#>.mat
  

%%%%%%%%% get ROIs for selected plane, view, and channel %%%%%%%%
% by using SVD file computed above. 
update_roi(ops0,handles)
% ROI file is stored in  
% ROI_<mouse_name>_<date>_plane<#>_ch<#>.mat
     

%%%%%%%%% get signals for selected plane, view, channel %%%%%%%%
% by using ROI file and registered Tiff computed above, 
update_signal(ops0,handles);
% Signal file is stored in
% Fsig_<mouse_name>_<date>_plane<#>_ch<#>.mat

%%%%% cleanup %%%%%
try
    for ii=1:length(ops1)
        if ops1{ii}.DeleteBin
            fclose('all');
            delete(ops1{ii}.RegFile);        % delete temporary bin file
        end
    end
end

function ops1=update_reg(ops0)

%%%%%%%%%%%% image registration and svd computation %%%%%%%%%%%%
ops1 = [];
if ops0.doRegistration
     ops1         = reg2P_kh004(ops0);  % do registration. 
     % At this point, ops1 contains Plane x View x Channel size. 
    % get_SVD function relies on temporally generated temp_file, so always use right after reg2P. 
     ops1         = get_svd_20171102(ops1); % get svd 
end

function ops1 = update_roi(ops0,handles)
%%%%%%%%% get ROIs for selected plane, view, and channel %%%%%%%%
PlaneChString = get(handles.lb_PlaneCh,'String');
PlaneChString =  PlaneChString(get(handles.lb_PlaneCh,'Value'));

if ops0.getROIs
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
       ops1{pp} = update_ops(ops0,tmp.ops);
       ops1{pp} = get_roi_20171102(ops1{pp},tmp.U,tmp.Sv,handles.DB.clustrules);
    end
    
end

function ops1 = update_signal(ops0,handles)
%%%%%%%%% get signals for selected plane, view, channel %%%%%%%%
% by using ROI file and registered Tiff computed above two. 

PlaneChString = get(handles.lb_PlaneCh,'String');
PlaneChString =  PlaneChString(get(handles.lb_PlaneCh,'Value'));
if ops0.DoGetSignal
       ops1 = [];
  for pp=1:length(PlaneChString)
      MinClustSize = handles.DB.clustrules.MinClust;
      ROIFile = sprintf('%s/ROI_%s_%s_%s_MinClust%d.mat', ops0.ResultsSavePath, ...
          ops0.mouse_name, ops0.date, PlaneChString{pp},MinClustSize);
     
      if exist(ROIFile,'file')
          fprintf('Loading ROI file: %s\n', ROIFile);
          tmp=load(ROIFile,'ops', 'res', 'stat', 'stat0', 'res0', 'clustrules');

          % tmp.ops contains previous options. ops0 inherits latest options.
          % ops1{pp} needs to be constructed based on ops0, but otherwise use old values
          ops1{pp} = update_ops(ops0,tmp.ops);
          
          ops1{pp}.PlaneChString = PlaneChString{pp};
          get_signal_20171102(ops1{pp});
          
      else
          errordlg(sprintf('%s does not exist.', ROIFile));
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
