function varargout = RoiGui_007(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RoiGui_007_OpeningFcn, ...
                   'gui_OutputFcn',  @RoiGui_007_OutputFcn, ...
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

%% data dependencies 
%  
%   h.dat.stat=update_stat(h.dat.res.iclust1,h.dat.res.M); 
%  h=init_cl(h);     % cluster selection and related variables
%    h=update_Ftrace(h); % inside, h.dat.cl.statTBL.skewF is updated.
%%
% --- Executes just before RoiGui_007 is made visible.
function RoiGui_007_OpeningFcn(hObject, eventdata, h, varargin)
% This function has no output args, see OutputFcn.
% varargin   command line arguments to RoiGui_007 (see VARARGIN)
 
h=init_icons(hObject,eventdata,h);
h=init_movie_timer(hObject,eventdata,h);
h.control_on = 0;
h.shift_on = 0;
 
h.output = hObject;
guidata(hObject, h);

% if nargin>=4 
%     pb_LoadPlane_Callback(handles.pb_loadPlane, eventdata, h,varargin);
% end

% UIWAIT makes RoiGui_007 wait for user response (see UIRESUME)
% uiwait(h.figure1);

function h=init_icons(hObject,eventdata,h)
h.icons.SVM_BW  =imread('SVM_Black_26x26.png');
h.icons.SVM_Color =imread('SVM_26x26.png');
h.icons.ROICircle_BW=imread('ROICircle_25x26_Gray.png');
h.icons.ROICircle_Color=imread('ROICircle_25x26_Color.png');
h.icons.CellDisplay_BW=imread('CellDisplay_25x26_Gray.png');
h.icons.CellDisplay_Color=imread('CellDisplay_25x26_Color.png');


function h=init_movie_timer(hObject,eventdata,h)
h.movie_timer = timer ;
h.movie_timer.Period = 0.05 ;
h.movie_timer.ExecutionMode = 'fixedSpacing' ;
h.movie_timer.TimerFcn = {@movie_timer_calback ,hObject};
% h.movie_timer.userdata = hObject;

function movie_timer_calback(obj,event,hObject)
h=guidata(hObject);

%  h=guidata(obj.userdata);
 
ind=get(h.axes_left,'UserData'); 
% first index is the current image index, second is the stride of the movie,  role forward (+1) or backward (-1).

nextind = ind(1)+ind(2); 
if nextind <1
%     nextind =  size(h.dat.reg_data,3);
    nextind =  h.dat.reg_data.TotalNSeries
end
if nextind >  h.dat.reg_data.TotalNSeries
    nextind = 1;
end
% 
LeftH = h.left_imageH;

img = h.dat.reg_data.get_single_slice(nextind);
img=img(h.dat.ops.yrange,h.dat.ops.xrange);
% img=repmat(img,[1,1,3]);

% img(:,:,1)=img(:,:,1)+h.dat.cellmask;b
img = imfuse(img,max(img(:))*uint16(h.dat.cellmask),'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
set(LeftH,'CData',img);

set(h.axes_left,'UserData',[nextind,ind(2)]);
set(h.uipanel_bottom,'Title',sprintf('(%d/%d)',nextind,h.dat.reg_data.TotalNSeries));
set(h.plot_fluorescence_timing_bar,'XData',nextind*h.dat.reg_decimation*[1 1],...
    'YData',get(h.axes_fluorescence,'YLim'));
% 
drawnow;
% disp('running')
    
% --- Outputs from this function are returned to the command line.
function varargout = RoiGui_007_OutputFcn(hObject, eventdata, h) 
% varargout  cell array for returning output args (see VARARGOUT);
varargout{1} = h.output;

function h=pb_LoadPlane_Callback(hObject, eventdata, h,varargin)
% [filename1,filepath1]=uigetfile('\\zserver\Lab\Share\Marius\', 'Select Data File');

if nargin>=4
    LOADNAME=varargin{1};
else
    LOADNAME=[];
end

if isfield(h, 'dat') && isfield(h.dat, 'filename')
    root = fileparts(h.dat.filename);
else
    root = 'G:Kosuk\DATA\F\';
end

if isempty(LOADNAME)
    [filename1,filepath1]=uigetfile(root, 'Select Data File');
    LOADNAME = fullfile(filepath1, filename1);
end

fprintf('Loading %s...\n',LOADNAME);
% W=whos('-file',LOADNAME);

h.dat = load(LOADNAME);
h.dat.filename = LOADNAME;


[filepath1,filename1,ext]=fileparts(LOADNAME);
set(h.figure1, 'Name', filename1);

%% initialize images 
% delete(h.axes_left.Children);
%% 


% if flag
    % if the user selected a file, do all the initializations
    rng('default')


% keyboard;
if isfield(h.dat, 'dat')
    h.dat = h.dat.dat;
    fprintf('h.dat.dat -> h.dat\n')
    h.dat.stat=update_stat(h.dat.res.iclust,h.dat.res.M); % update stat at loading point.
%     h=update_Ftrace(h);
%     h=update_cl(h);     % cluster selection and related variables
   
elseif isfield(h.dat.ops,'processed_date')
    fprintf('processed data \n')
    if ~isfield(h.dat.res,'iclust1')
     h.dat.res.iclust1 = h.dat.res.iclust;
    end
 
    h.dat.stat=update_stat(h.dat.res.iclust,h.dat.res.M); % update stat at loading point.
      h=init_QV(h);
%     h=update_Ftrace(h);
%     h=update_cl(h);     % cluster selection and related variables
    
else % first time to load proc. 
    
    h.dat.filename = getOr(h.dat,'filename',fullfile(filepath1, filename1));
    
%     if isfield(h.dat.stat, 'parent')
%         h = get_parent_stats(h);
%     end
    
    %h.dat.res.iclust = reshape(h.dat.res.iclust, h.dat.res.Ly, h.dat.res.Lx);
    h.dat.res.iclust = getOr(h.dat.res,'iclust',...
                       reshape(h.dat.res.iclust, h.dat.res.Ly, h.dat.res.Lx));
    h.dat.res.iclust1 = h.dat.res.iclust; % update iclust to merge. 
    

    h.dat.stat=update_stat(h.dat.res.iclust1,h.dat.res.M); % stat: statistics or features of each cluster (ROI)
    h.dat.ops.Nk = numel(h.dat.stat);
    h=init_cl(h);     % cluster selection and related variables
   
    if isfield(h.dat, 'clustrules')
         % ROI rules
         if ~isfield(h.dat.clustrules,'MinNpix')
             h.dat.clustrules=get_clustrules(h.dat.clustrules);
         end
         h.dat.res.Mrs_thresh_orig = h.dat.clustrules.Compact;
         h.dat.cl.npix_low_orig    = h.dat.clustrules.MinNpix;
         h.dat.cl.npix_high_orig   = h.dat.clustrules.MaxNpix;
    else
        % ROI rules
        h.dat.res.Mrs_thresh_orig   = 3;
        h.dat.cl.npix_low_orig      = 20;
        h.dat.cl.npix_high_orig     = 500;
    end
       
    h = setOriginalThresh(h);
    
    set(h.edit_Compactness,'String', num2str(h.dat.res.Mrs_thresh));
    set(h.edit_PixelCountHigh,'String', num2str(h.dat.cl.npix_high));
    set(h.edit_PixelCountLow,'String', num2str(h.dat.cl.npix_low));
    
    % quadrant view 
    h=init_QV(h);
    
    % start with unit vector map
    h.dat.display_select = true; % start with selection view
    
     
    Mmin = min(h.dat.res.M(:));
    Mmax = max(h.dat.res.M(:));
    h.dat.res.M0 = (h.dat.res.M-Mmin+0.05)/(Mmax-Mmin);
        
    % <To Do>: check the meaning of this
    % loop through redcells and set h.dat.cl.rands(h.dat.F.ichosen) = 0
%     for j = find(h.dat.cl.redcell)
%         if ~isempty(j)
%             h.dat.F.ichosen = j;
%             h.dat.cl.rands(h.dat.F.ichosen) = 0;
%         end
%     end
    
    icell = find(h.dat.cl.selected);
    if ~isempty(icell)
        h.dat.F.ichosen = icell(1); %ceil(rand * numel(icell))
    else
        h.dat.F.ichosen = 1; %ceil(rand * numel(icell))
    end
    
    % Fcell is raw F within ROI
    h.dat.F.Fcell = h.dat.Fcell; 
    h.dat=rmfield(h.dat,'Fcell');   
    
    % FcellNeu is raw F within Neuropil Mask. 
    h.dat.F.FcellNeu = h.dat.FcellNeu; 
    h.dat = rmfield(h.dat,'FcellNeu');
    
    for cc=1:length(h.dat.F.Fcell)
        Coef='1/Regress';
        h.dat.F.Ftrue{cc} =subtract_neurop(h.dat.F.Fcell{cc},...
            h.dat.F.FcellNeu{cc},Coef);
    end
end

%% added 20210118 to make sure all the scalar stats are scalars in vector, not in cells. -> should not be necessary if statTBL is correctly updated after splitting and merging in get_stat_from_iclust
% scalar_vars_in_statTBL={'mrs','npix','mrs0','Compactness', 'V','Solidity','Eccentricity','Perimeter','SNratio','skewF'};
% for vv=1:length(scalar_vars_in_statTBL)
%     varname=scalar_vars_in_statTBL{vv};
%     if iscell(h.dat.cl.statTBL.(varname))
%         varval=[h.dat.cl.statTBL.(varname){:}]; % we found some statTBL values are empty, and that's why it became cell array to keep it empty.
%         h.dat.cl.statTBL.(varname)=varval(:);
%     end
% end
%%
h.dat.graph.maxmap = 1;
ops = h.dat.ops;

% h.dat.graph.mimg(:,:,1)  is mean image
if isfield(ops, 'mimg1') && ~isempty(ops.mimg1)
    h.dat.graph.maxmap = h.dat.graph.maxmap + 1;
    h.dat.graph.mimg(:,:,h.dat.graph.maxmap) = ops.mimg1(ops.yrange, ops.xrange);
    h.dat.graph.mimg_proc(:,:,h.dat.graph.maxmap) = normalize_image(h.dat.graph.mimg(:,:,h.dat.graph.maxmap));
end

if isfield(ops, 'mimgRED') && ~isempty(ops.mimgRED)
    h.dat.graph.maxmap = h.dat.graph.maxmap + 1;
    h.dat.graph.mimg(:,:,h.dat.graph.maxmap) = ops.mimgRED(ops.yrange, ops.xrange);
    h.dat.graph.mimg_proc(:,:,h.dat.graph.maxmap) = normalize_image(h.dat.graph.mimg(:,:,h.dat.graph.maxmap));
end

if isfield(ops, 'mimgREDcorrected') && ~isempty(ops.mimgREDcorrected)
    h.dat.graph.maxmap = h.dat.graph.maxmap + 1;
    h.dat.graph.mimg(:,:,h.dat.graph.maxmap) = ops.mimgREDcorrected;
    h.dat.graph.mimg_proc(:,:,h.dat.graph.maxmap) = normalize_image(h.dat.graph.mimg(:,:,h.dat.graph.maxmap));
end

h.dat.procmap = 0;

h=update_cl(h);    
h=update_Ftrace(h); % inside, h.dat.cl.statTBL.skewF is updated.
h.dat.stat=table2struct(h.dat.cl.statTBL);

h = ApplyROIFilter(h);
h.dat.F.ichosen_append = [];
h = buildSat(h);
h = buildHue(h);
h = buildLambdaValue(h);
    
% to init image in left/right axes
axes(h.axes_left); 
h.left_imageH=imshow(h.dat.graph.mimg(:,:,h.dat.graph.map)); 
colormap('gray'); axis off
% set(h.left_imageH,'CDataMapping','direct');

% axes(h.axes_right);  
% h.left_imageH=imagesc(h.dat.graph.mimg(:,:,h.dat.graph.map)); 
% colormap('gray'); axis off
axes(h.axes_fluorescence); 

[NN NT] = size(h.dat.F.trace);

% in case there remains previous drawings
 try 
     delete(get(h.axes_fluorescence,'Children'))
 end

 h.plot_fluorescence_F = init_fluorescence_plot_handles(h, 100); % initially, prepare 100 lines
  
% h.plot_fluorescence_F = plot(1:NT,zeros(1,NT));hold on;
% h.plot_fluorescence_F_blue = plot(1:NT,zeros(1,NT),'b');
% h.plot_fluorescence_F_red = plot(1:NT,zeros(1,NT),'r');

h.plot_fluorescence_baseF = plot(1:NT,zeros(1,NT));
h.plot_fluorescence_FCell = plot(1:NT,zeros(1,NT));
h.plot_fluorescence_timing_bar = line([0 0], [-1 1]);

% init Display mode. 
set(h.tg_DisplayMode,'Value',1);


h= ModeSelectionButtonGroup_SelectionChangedFcn(h.tg_DisplayMode, eventdata, h);
 
ButtonDownFcn=@(hObject,eventdata)RoiGui_007('TimeSelectionButtonDownFcn',h.axes_fluorescence,eventdata,h);
set(h.axes_fluorescence,'ButtonDownFcn',ButtonDownFcn);

h=renew_cluster(h);
h=redraw_fluorescence(h);
h=redraw_figure(h);

guidata(hObject,h)

% end

function  stat=update_stat(iclust,M)
% L = h.dat.res.iclust;
% M = h.dat.res.M;
fprintf('Updating stat...\n');
stat =  get_stat_from_iclust(iclust,M);

function h=init_cl(h)

for ii=1:length(h.dat.stat)
    h.dat.stat(ii).Compactness=h.dat.stat(ii).mrs/h.dat.stat(ii).mrs0;
end

N=length(h.dat.stat);
h.dat.cl.statTBL = struct2table(h.dat.stat);
h.dat.cl.selected = ones(N,1);
h.dat.cl.statTBL.skewF = inf(N,1);
h.dat.cl.statTBL.SNratio = ones(N,1); 

rng('default');
h.dat.cl.rands_orig   = .1 + .8 * rand(h.dat.ops.Nk,1);
h.dat.cl.rands        = h.dat.cl.rands_orig;
h.dat.cl.manualmerge  = getOr(h.dat.cl,'manualmerge',zeros(h.dat.ops.Nk,1));
h.dat.cl.svm_class    = getOr(h.dat.cl,'svm_class',zeros(h.dat.ops.Nk,1));
h.dat.cl.svm_predicted_type   = zeros(h.dat.ops.Nk,1);
h.dat.cl.manual  = getOr(h.dat.cl,'manual',zeros(h.dat.ops.Nk, 1));
h.dat.cl.redcell = getOr(h.dat.cl,'redcell',zeros(h.dat.ops.Nk, 1));
 

% roi_type, 'type', RGB
h.dat.cl.type_number_table = ...
    {0,'Unknown',   [1 1 1];...
    1, 'Cell'       [1 0 0]; ...
    2,'Dendrite',   [0 1 0];...
    3, 'Axon',      [0 0 1];...
    4, 'Spine',     [1 1 0];...
    9, 'Noise',     [1 1 1]};
    
function h=update_cl(h)

% should be no longer necessary, as it is taken care by get_stat_from_iclust 

% for ii=1:length(h.dat.stat)
%     if isempty(h.dat.stat(ii).Compactness)
%         h.dat.stat(ii).Compactness=h.dat.stat(ii).mrs/h.dat.stat(ii).mrs0;
%     end
% end

h.dat.cl.statTBL = struct2table(h.dat.stat);
N = length(h.dat.stat);
% number of new entries, cl.selected is not updated yet. 
dN = N-length(h.dat.cl.selected);

if dN>0
    
VarsToUpdate   = setdiff(h.dat.cl.statTBL.Properties.VariableNames,fieldnames(h.dat.stat));

if ~isempty(VarsToUpdate)
    % rescue variables that are not included in h.dat.stat but in h.dat.cl.statTBL
    for ss=VarsToUpdate'
        tmp.(ss{1})= h.dat.cl.statTBL.(ss{1});
    end
end

h.dat.cl.selected      = cat(1,h.dat.cl.selected,ones(dN,1));
    if ~isempty(VarsToUpdate)
 
    for ss=VarsToUpdate'
        switch ss{1}
            case 'skewF'
                h.dat.cl.statTBL.(ss{1})= cat(1,tmp.(ss{1}),inf(dN,1));
            case 'SNratio'
                h.dat.cl.statTBL.(ss{1})= cat(1,tmp.(ss{1}),ones(dN,1));
            otherwise
                error('Unknown variables %s',ss);
        end
    end
    end
    
    rng('default');
    h.dat.cl.rands_orig   = .1 + .8 * rand(N,1);
    h.dat.cl.rands        = h.dat.cl.rands_orig;
    h.dat.cl.manualmerge  = cat(1,h.dat.cl.manualmerge,zeros(dN,1));
    h.dat.cl.svm_class    = cat(1,h.dat.cl.svm_class,  zeros(dN,1));
    h.dat.cl.svm_predicted_type   = cat(1,h.dat.cl.svm_predicted_type,  zeros(dN,1));
    h.dat.cl.manual       = cat(1,h.dat.cl.manual,ones(dN,1));
    h.dat.cl.redcell      = cat(1,h.dat.cl.redcell,zeros(dN,1));

end

function h=add_cl(h,dN)
% add new rows of statTBL for split or merged ROI. 

% should be no longer necessary, as it is taken care by get_stat_from_iclust 

% for ii=1:length(h.dat.stat)
%     if isempty(h.dat.stat(ii).Compactness)
%         h.dat.stat(ii).Compactness=h.dat.stat(ii).mrs/h.dat.stat(ii).mrs0;
%     end
% end

h.dat.cl.statTBL = struct2table(h.dat.stat);

if dN>0
    
VarsToUpdate   = setdiff(h.dat.cl.statTBL.Properties.VariableNames,fieldnames(h.dat.stat));

if ~isempty(VarsToUpdate)
    % rescue variables that are not included in h.dat.stat but in h.dat.cl.statTBL
    for ss=VarsToUpdate'
        tmp.(ss{1})= h.dat.cl.statTBL.(ss{1});
    end
end

h.dat.cl.selected      = cat(1,h.dat.cl.selected,ones(dN,1));
%     if ~isempty(VarsToUpdate)
 
    for ss=VarsToUpdate'
        switch ss{1}
            case 'skewF'
                h.dat.cl.statTBL.(ss{1})= cat(1,tmp.(ss{1}),inf(dN,1));
            case 'SNratio'
                h.dat.cl.statTBL.(ss{1})= cat(1,tmp.(ss{1}),ones(dN,1));
            otherwise
                error('Unknown variables %s',ss);
        end
    end
%     end
    
    rng('default');
    h.dat.cl.rands_orig   = .1 + .8 * rand(N,1);
    h.dat.cl.rands        = h.dat.cl.rands_orig;
    h.dat.cl.manualmerge  = cat(1,h.dat.cl.manualmerge,zeros(dN,1));
    h.dat.cl.svm_class    = cat(1,h.dat.cl.svm_class,  zeros(dN,1));
    h.dat.cl.svm_predicted_type   = cat(1,h.dat.cl.svm_predicted_type,  zeros(dN,1));
    h.dat.cl.manual       = cat(1,h.dat.cl.manual,ones(dN,1));
    h.dat.cl.redcell      = cat(1,h.dat.cl.redcell,zeros(dN,1));

end

function h=init_QV(h)
 % set all quadrants as not visited
    h.QV.quadvalue = zeros(3);
    for j = 1:3
        for i = 1:3
            set(h.(sprintf('Q%d%d', j,i)), 'BackgroundColor',[.92 .92 .92]);
        end
    end
    h.dat.QV.ylim = [0 h.dat.res.Ly];
    h.dat.QV.xlim = [0 h.dat.res.Lx];    
    % x and y limits on subquadrants
    h.dat.QV.x0all = round(linspace(0, 19/20*h.dat.res.Lx, 4));
    h.dat.QV.y0all = round(linspace(0, 19/20*h.dat.res.Ly, 4));
    h.dat.QV.x1all = round(linspace(1/20 * h.dat.res.Lx, h.dat.res.Lx, 4));
    h.dat.QV.y1all = round(linspace(1/20 * h.dat.res.Ly, h.dat.res.Ly, 4));
    

function obj=init_fluorescence_plot_handles(h, n)
obj = zeros(n,1);

[NN NT] = size(h.dat.F.trace);
for ii=1:n
    obj(ii) = plot(1:NT,zeros(1,NT));
    hold on;
end


function h=update_Ftrace(h)
%% construct F.trace
h.dat.graph.map = 1;
h.dat.F.trace = [];
h.dat.F.truetrace = [];
for i = 1:length(h.dat.F.Fcell)
    h.dat.F.trace = cat(2, h.dat.F.trace, h.dat.F.Fcell{i});
%     h.dat.F.truetrace = cat(2, h.dat.F.truetrace, h.dat.F.Ftrue{i});
end

if isfield(h.dat.F, 'FcellNeu')
    h.dat.F.neurop = [];
    for i = 1:length(h.dat.F.FcellNeu)
        h.dat.F.neurop = cat(2, h.dat.F.neurop, h.dat.F.FcellNeu{i});
    end    
    
else
   h.dat.F.neurop = zeros(size(h.dat.F.trace), 'single');
end

 Coef='1/Regress';
h.dat.F.truetrace=subtract_neurop(h.dat.F.trace,h.dat.F.neurop,Coef);
 
    
h.dat.plot_neu = 0;
% Is skewF calculated from baseline subtracted data? or raw data? 

use_trueF_for_skewF = 1;
if (use_trueF_for_skewF)
    % true F version
    [skewF, SNratio,C_of_zF, C_of_dzF]=get_F_dependent_stats(h.dat.F.truetrace);
else
    % raw trace version
    [skewF, SNratio,C_of_zF, C_of_dzF]=get_F_dependent_stats(h.dat.F.trace);
end

dN=height(h.dat.cl.statTBL)-length(skewF);
if dN>0
    skewF=cat(1,skewF(:),inf(dN,1));
    SNratio=cat(1,SNratio(:),ones(dN,1));
end
h.dat.cl.statTBL.skewF =skewF;
h.dat.cl.statTBL.SNratio  = SNratio;

h.dat.cl.C_of_zF=C_of_zF;
h.dat.cl.C_of_dzF=C_of_dzF;
h.dat.cl.IDX = recalc_IDX(h);


% if (use_trueF_for_skewF)
%     % true F version
%     h.dat.cl.statTBL.skewF = skewness(h.dat.F.trace,0,2);
% %     h.dat.cl.statTBL.AutoCorr = fft(h.dat.F.trace,[],2);
% else
%     % raw trace version
%     h.dat.cl.statTBL.skewF = skewness(h.dat.F.truetrace,0,2);
% %     h.dat.cl.statTBL.AutoCorr = fft(h.dat.F.truetrace,[],2);
% end
% 
% zF = zscore(h.dat.F.truetrace,0,2);
% d = 10;
% dzF = zscore(zF(:,(1+d):end)-zF(:,1:end-d),0,2);
% % cellid = find(h.dat.cl.selected); % to visualize correlation, 
% % h.dat.cl.C_of_zF=zF(cellind,:)*zF(cellind,:)'/size(zF,2);
% % would be easier to see.  
% h.dat.cl.C_of_zF=zF*zF'/size(zF,2);
% h.dat.cl.C_of_dzF=dzF*dzF'/size(dzF,2);
% h.dat.cl.IDX = recalc_IDX(h);
% 
% %% Here is a new code to add signal/noise ratio by calculating the power ratio below 1Hz and above. 
% FPS= 30;
% Tlen = 100*FPS;
% % Multitaper Time-Frequency Power-Spectrum (power spectrogram)
% % function A=mtpsg(x,nFFT,Fs,WinLength,nOverlap,NW,nTapers)
% % x : input time series
% nFFT = 2^nextpow2(Tlen); %number of points of FFT to calculate (default 1024)
% Fs = FPS; %sampling frequency (default 2)
% WinLength = nFFT; %length of moving window (default is nFFT)
% nOverlap = 0;%nFFT/2; %overlap between successive windows (default is WinLength/2)
% NW = 3; %time bandwidth parameter (e.g. 3   or 4), default 3
% nTapers = 2*NW-1; % nTapers = number of data tapers kept, default 2*NW -1
% Detrend= 1; 
% 
% SNratio = nan(1,size(zF,1));
% 
% waitH = waitbar(0,'Computing power spectrum ...');
% for ii=1:size(zF,1)
% [yo, fo]=my_mtcsg(zF(ii,:),nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers);
% myo=mean(abs(yo),2);
% % plot(fo,myo);
% 
% ind=fo<1; 
% sig1Hz_power=sum(myo(ind));
% noise_power=sum(myo(~ind));
% SNratio(ii)=sig1Hz_power/noise_power;
% waitbar(ii/size(zF,1),waitH);
% end
% close(waitH);
% SNratio(isnan(SNratio))=1;
% h.dat.cl.statTBL.SNratio = log(SNratio(:));


function [cell_index,roi_index]=get_correlated_roi(C_of_zF,cell_id,threshold)
%  [cell_index,roi_index]=get_correlated_roi(C_of_zF,cell_id,threshold);
%  Input: 
%  C_of_zF:  cross-covariance or normalized cross-correlation of F. The value is between -1 to 1. 
%  cell_id:  cell id within iscell=1 group. 
%  threshold: Correlation higher than this value is returned.

cell_index=find(C_of_dzF(cell_id,:)>threshold);




% function pushbutton61_Callback(hObject, eventdata, h)
% % keep TOP variance region
% h.dat.cl.topregion = zeros(h.dat.res.Ly, h.dat.res.Lx);
% xs = repmat(1:h.dat.res.Lx, h.dat.res.Ly, 1);
% ys = repmat((1:h.dat.res.Ly)', 1, h.dat.res.Lx);
% 
% for k = 1:length(h.dat.stat)
%     if ~isempty(h.dat.stat(k).Vregion)
%         [~, itop] = max(h.dat.stat(k).Vregion);
%         h.dat.cl.topregion(h.dat.stat(k).region{itop}) = 1;
%         h.dat.cl.statTBL.npix(k) = h.dat.stat(k).npixels(itop);
%     end
% end
% 
% h = buildLambdaValue(h);
% h=redraw_figure(h);
% guidata(hObject,h)

% function pb_binarymask_Callback(hObject, eventdata, h)
% % binary mask
% h.dat.graph.img0.V = ones(h.dat.res.Ly, h.dat.res.Lx);
% h = buildLambdaValue(h);
% 
% guidata(hObject,h);
% h=redraw_figure(h);


% function pb_variancemask_Callback(hObject, eventdata, h)
% % variance explained mask
% h.dat.graph.img0.V = reshape(h.dat.res.M, h.dat.res.Ly, h.dat.res.Lx)/h.dat.cl.MeanM;
% h = buildLambdaValue(h);
% guidata(hObject,h);
% h=redraw_figure(h);


function pb_defaulthue_Callback(hObject, eventdata, h)
% original default hue
h.dat.cl.rands   = h.dat.cl.rands_orig;
h.dat.graph.img1.H       = reshape(h.dat.cl.rands(h.dat.res.iclust1), h.dat.res.Ly, h.dat.res.Lx);
h.dat.graph.img2.H       = reshape(h.dat.cl.rands(h.dat.res.iclust1), h.dat.res.Ly, h.dat.res.Lx);
guidata(hObject,h);
h=redraw_figure(h);

function pb_randomizehue_Callback(hObject, eventdata, h)
% randomize hue
rng('shuffle') 
unique_iclust = unique(h.dat.res.iclust1(:));
h.dat.cl.rands     = rand(1, max(unique_iclust));
h.dat.cl.rands(1)  = .15;
h.dat.graph.img1.H       = reshape(h.dat.cl.rands(h.dat.res.iclust1), h.dat.res.Ly, h.dat.res.Lx);
h.dat.graph.img2.H       = reshape(h.dat.cl.rands(h.dat.res.iclust1), h.dat.res.Ly, h.dat.res.Lx);

guidata(hObject,h);
h=redraw_figure(h);



% function edit33_Callback(hObject, eventdata, h)
% h.dat.cl.pixthresh_percent = str2double(get(h.edit33,'String'));
% % h = exclude_pixels_percent(h);
% h = buildLambdaValue(h);
% h=redraw_figure(h);
% guidata(hObject,h);

% function h = exclude_pixels_percent(h)
% h.dat.cl.excl_pix_perc = zeros(h.dat.res.Ly, h.dat.res.Lx);
% for k = 1:h.dat.ops.Nk
%     which_pix = find(h.dat.res.iclust==k);
%    [Msort, isort] = sort(h.dat.res.M(which_pix), 'ascend'); 
%    Msort = cumsum(Msort);
%    Msort = 100 * Msort/max(Msort);
%    ifi = find(Msort>h.dat.cl.pixthresh_percent, 1);
%    h.dat.cl.excl_pix_perc(which_pix(isort(1:ifi))) = 1;
% end


function edit_Compactness_Callback(hObject, eventdata, h)
h.dat.res.Mrs_thresh = str2double(get(h.edit_Compactness,'String'));
h = ApplyROIFilter(h);
h = buildLambdaValue(h);
h=redraw_figure(h);
guidata(hObject,h);

function edit_Compactness_CreateFcn(hObject, eventdata, h)
set(hObject,'String', num2str(2));
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function pb_saveprocfile_Callback(hObject, eventdata, h)

if ~isempty(strfind(h.dat.filename,'Spk'))
    filename = h.dat.filename;
else
    [filepath,filename,ext]=fileparts(h.dat.filename);
    filename = strrep(filename,'_proc','');
    filename = fullfile(filepath,[filename,'_proc.mat']);
end
fprintf('Saving results %s ...\n',filename)


dat = h.dat;
dat.res.iclust = dat.res.iclust1; % overwrite iclust 
% this eliminates some split ROI's area to be zero, and makes it difficult to calculate the later process. 
 dat.stat=update_stat(dat.res.iclust,dat.res.M); % <We should update it at loading point>

dat.F.trace = []; % no longer needed
dat.F.truetrace =[]; % no longer needed
dat.F.neurop =[]; % no longer needed
dat.F.Ftrue =[]; % no longer needed
% dat.res = rmfield(dat.res,'iclust1');

if isfield(dat,'reg_data'),dat=rmfield(dat,'reg_data'); end % to remove huge data. end
if isfield(dat,'svd'),dat=rmfield(dat,'svd'); end % to remove huge data. end

dat.ops.processed_date = datestr(now,'yyyy_mm_dd_HH_MM_SS');
dat.ops.processed_mfile = mfilename; 
dat.ops.ProcFileName = filename; 

% try
    save(filename, '-struct','dat')
% catch
%     save(filename,'-v7.3', 'dat')
% end
         
            
pb_Print_Callback(hObject, eventdata, h);

% print(gcf,'-dpdf','-bestfit',[printname,'.pdf']);

function edit_PixelCountHigh_Callback(hObject, eventdata, h)
h.dat.cl.npix_high = str2double(get(h.edit_PixelCountHigh,'String'));
h = ApplyROIFilter(h);
h = buildLambdaValue(h);
h=redraw_figure(h);
guidata(hObject,h);

function edit_PixelCountHigh_CreateFcn(hObject, eventdata, h)
set(hObject,'String', num2str(400));
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_PixelCountLow_Callback(hObject, eventdata, h)
h.dat.cl.npix_low = str2double(get(h.edit_PixelCountLow,'String'));
h = ApplyROIFilter(h);
h = buildLambdaValue(h);
h=redraw_figure(h);
guidata(hObject,h);

function edit_PixelCountLow_CreateFcn(hObject, eventdata, h)
set(hObject,'String', num2str(30));
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function figure1_ResizeFcn(hObject, eventdata, h)


% --- Executes on button press in Q11.
function Q11_Callback(hObject, eventdata, h)
iy = 1; ix = 1;
quadrant(hObject, h, iy, ix)
paint_quadbutton(h, iy, ix);

% --- Executes on button press in Q12.
function Q12_Callback(hObject, eventdata, h)
iy = 1; ix = 2;
quadrant(hObject, h, iy, ix)
paint_quadbutton(h, iy, ix);

% --- Executes on button press in Q13.
function Q13_Callback(hObject, eventdata, h)
iy = 1; ix = 3;
quadrant(hObject, h, iy, ix)
paint_quadbutton(h, iy, ix);

% --- Executes on button press in Q21.
function Q21_Callback(hObject, eventdata, h)
iy = 2; ix = 1;
quadrant(hObject, h, iy, ix)
paint_quadbutton(h, iy, ix);

% --- Executes on button press in Q22.
function Q22_Callback(hObject, eventdata, h)
iy = 2; ix = 2;
quadrant(hObject, h, iy, ix)
paint_quadbutton(h, iy, ix);

% --- Executes on button press in Q23.
function Q23_Callback(hObject, eventdata, h)
iy = 2; ix = 3;
quadrant(hObject, h, iy, ix)
paint_quadbutton(h, iy, ix);

% --- Executes on button press in Q31.
function Q31_Callback(hObject, eventdata, h)
iy = 3; ix = 1;
quadrant(hObject, h, iy, ix)
paint_quadbutton(h, iy, ix);

% --- Executes on button press in Q32.
function Q32_Callback(hObject, eventdata, h)
iy = 3; ix = 2;
quadrant(hObject, h, iy, ix)
paint_quadbutton(h, iy, ix);

% --- Executes on button press in Q33.
function Q33_Callback(hObject, eventdata, h)
iy = 3; ix = 3;
quadrant(hObject, h, iy, ix);
paint_quadbutton(h, iy, ix);

function paint_quadbutton(h, iy, ix)
   
    for j = 1:3
        for i = 1:3
            if h.QV.quadvalue(j,i)==1
                set(h.(sprintf('Q%d%d', j,i)), 'BackgroundColor','yellow');
            end
        end
    end
    
    set(h.(sprintf('Q%d%d', iy,ix)), 'BackgroundColor','red');
    
% --- Executes on button press in full.
function full_Callback(hObject, eventdata, h)
%QV: quadrant view
h.dat.QV.ylim = [0 h.dat.res.Ly];
h.dat.QV.xlim = [0 h.dat.res.Lx];
guidata(hObject,h);
h=redraw_figure(h);

function quadrant(hObject, h, iy, ix)
h.dat.QV.ylim = [h.dat.QV.y0all(iy) h.dat.QV.y1all(iy+1)];
h.dat.QV.xlim = [h.dat.QV.x0all(ix) h.dat.QV.x1all(ix+1)];
h.QV.quadvalue(iy, ix) = 1;

guidata(hObject,h);
h=redraw_figure(h);

% --- Executes on key press with focus on figure1 or any of its controls.
function figure1_WindowKeyPressFcn(hObject, eventdata, h)


switch eventdata.Key
      case 'control'
          h.control_on = 1;
          
%           fprintf('%s pressed.\n',eventdata.Key);
%     otherwise
%         error('unknown release')
end


if ~isempty(eventdata.Modifier)
    switch eventdata.Modifier{1}
        case {'shift','alt'}
            h.shift_on = 1;
            data.shift_on=h.shift_on;
            set(h.tg_DisplayMode,'UserData',data);
        case 'control'
            h.control_on = 1;
        otherwise
            h.shift_on = 0;
    end
end

h=ImageDisplaySwitchKeyPressFcn(hObject,eventdata,h);
guidata(hObject,h);


function h=ImageDisplaySwitchKeyPressFcn(hObject,eventdata,h)
% eventdata.Key
% eventdata.Modifier

if ~isempty(eventdata.Modifier)
    switch eventdata.Modifier{1}
        case {'shift','alt'}
            shift_on = 1;
            data.shift_on=shift_on;
            set(h.tg_DisplayMode,'UserData',data);
            
        case 'control'
            control_on = 0;
        otherwise
            shift_on = 0;
    end
end

switch eventdata.Key
     case 'a'
         if control_on
             fprintf('Select all\n');
         end
         
         
         if isempty(h.dat.F.ichosen_append), h.dat.F.ichosen_append=[];end
         
         Xrange=xlim;
         Yrange=ylim;
         Xrange=max(1,ceil(Xrange(1))):min(h.dat.res.Lx,floor(Xrange(end)));
         Yrange=max(1,ceil(Yrange(1))):min(h.dat.res.Ly,floor(Yrange(end)));
         
         
         
         if isempty(h.dat.F.ichosen_append), h.dat.F.ichosen_append=[];end

         selected_within_view=unique(h.dat.res.iclust1(Yrange,Xrange));
         selected_within_view(selected_within_view==1)=[]; % remove noise;
         
         h.dat.F.ichosen = intersect(find(h.dat.cl.selected),selected_within_view);
         h.dat.F.ichosen_append=h.dat.F.ichosen;
         
         
         h=update_figure(h);
         guidata(hObject,h);
    case 'f'
%         % flip currently selected unit
%         h.dat.cl.manual(h.dat.F.ichosen) = .5 - h.dat.cl.selected(h.dat.F.ichosen);
%         h = ApplyROIFilter(h);
%         h = buildLambdaValue(h);
%         guidata(hObject,h);
%         if h.dat.graph.maxmap==1
%             h=redraw_figure(h);
%         end
%          last_keypress = eventdata.Key;
    case 's'
        % manual selection of units is removed, instead, flip display mode
        % of selected or non-selected.
%         fprintf('display selection=%d\n',h.dat.display_select);
        h.dat.display_select = ~h.dat.display_select;
        if h.dat.display_select
            set(h.uipanel_Left,'Title','Selected ROIs');
        else
            set(h.uipanel_Left,'Title','Unselected ROIs');
        end
        
        h=pb_rois_display_Callback(hObject, eventdata, h);
%          fprintf('display selection=%d\n',h.dat.display_select);
        last_keypress = eventdata.Key;
    case 'q'
        if h.dat.display_select
            set(h.uipanel_Left,'Title','Selected ROIs');
        else
            set(h.uipanel_Left,'Title','Unselected ROIs');
        end
        h.dat.graph.map = 1;
        h=pb_rois_display_Callback(hObject, eventdata, h);
        last_keypress = eventdata.Key;
    case 'c'
         set(h.uipanel_Left,'Title','Correlated ROIs');
        h.dat.graph.map = 2;
        h=pb_CorrROI_display_Callback(hObject, eventdata, h);
        last_keypress = eventdata.Key;
        
    case 'w'
        set(h.uipanel_Left,'Title','weighted image of corr-pixel ("p" to proc on/off)');
        h.dat.graph.map = 2;
        h=pb_meandisplay_Callback(hObject, eventdata, h);
        last_keypress = eventdata.Key;
  
    case 'e'
       
        h.dat.graph.map = 3;
        if h.dat.graph.maxmap>2
            pb_Rede_display_Callback(hObject, eventdata, h);
        end
        last_keypress = eventdata.Key;
    case 'r'
        h.dat.graph.map = 3;
        if h.dat.graph.maxmap>3
            pb_Rede_display_Callback(hObject, eventdata, h);
        end
        last_keypress = eventdata.Key;
    case 'p'
        h=pb_proccdisplay_Callback(hObject, eventdata, h);
        last_keypress = eventdata.Key;
    case 'm'
        h=redraw_movie(hObject,eventdata,h);
        last_keypress = eventdata.Key;
    case 'leftarrow'
        previous_frame(hObject,eventdata,h);
        last_keypress = eventdata.Key;
    case 'rightarrow'
        next_frame(hObject,eventdata,h);
        last_keypress = eventdata.Key;
    case 'period' % shift + '.' = >, toggle ON/OFF movie forward.
        if shift_on
            userdata=get(h.axes_left,'UserData');
            if length(userdata)==1
                userdata(2)=1; % init 
            end
               previous_stride = userdata(2);
            if strcmp(h.movie_timer.Running,'on') && previous_stride >0  % if previously running forward and reaches max speed
                stop(h.movie_timer); % stop it.
            elseif strcmp(h.movie_timer.Running,'on') && previous_stride < 0 % if previously running backward
                userdata(2)=1; % role forward.
                set(h.axes_left,'UserData',userdata); % just change the stride.
            elseif strcmp(h.movie_timer.Running,'off') % if not started, start in forward direction.
                userdata(2)=1; % role forward.
                set(h.axes_left,'UserData',userdata);
                start(h.movie_timer)
            end
        end
        
        last_keypress= eventdata.Key;
    case 'comma' % shift + ',' = <, toggle ON/OFF movie backward.
        if shift_on
            userdata=get(h.axes_left,'UserData');
            previous_stride = userdata(2);
            if strcmp(h.movie_timer.Running,'on') && previous_stride <0 % if previously running backward
                stop(h.movie_timer); % stop it.
            elseif strcmp(h.movie_timer.Running,'on') && previous_stride >0 % if previously running forward
                userdata(2)=-1; % role backward.
                set(h.axes_left,'UserData',userdata); % just change the stride.
            elseif strcmp(h.movie_timer.Running,'off') % if not started, start in forward direction.
                userdata(2)=-1; % role backward.
                set(h.axes_left,'UserData',userdata);
                start(h.movie_timer)
            end
        end
        last_keypress= eventdata.Key;
        
    case 'slash' % shift + ',' = <, toggle ON/OFF movie backward.
        if shift_on
            userdata=get(h.axes_left,'UserData');
%             previous_stride = userdata(2);
            if strcmp(h.movie_timer.Running,'on')
                userdata(2)=userdata(2)*2; % double speed .
                set(h.axes_left,'UserData',userdata);
            end
        end
        last_keypress= eventdata.Key;
    otherwise
        last_keypress = [];
end

if ~isempty(last_keypress)
    h.last_keypress = last_keypress;
end
% h.dat.display_select
% guidata(hObject,h);

% --- Executes on key release with focus on figure1 or any of its controls.
function figure1_WindowKeyReleaseFcn(hObject, eventdata, h)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was released, in lower case
%	Character: character interpretation of the key(s) that was released
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) released
% handles    structure with handles and user data (see GUIDATA)
ImageDisplaySwitchKeyReleaseFcn(hObject,eventdata,h)

function ImageDisplaySwitchKeyReleaseFcn(hObject,eventdata,h)
%  eventdata.Key
%  eventdata.Modifier

%    eventdata
   
switch eventdata.Key
    case {'leftarrow', 'rightarrow'}
        if strcmp(h.movie_timer.Running,'on')
            stop(h.movie_timer)
        end
%     case 'control'
%         h.control_on = 0;
    case {'shift'}
     h.shift_on = 0;
%       fprintf('%s released.\n',eventdata.Key);
    case 'control'
          h.control_on = 0;
%           fprintf('%s released.\n',eventdata.Key);
%     otherwise
%         error('unknown release')
end


if ~isempty(eventdata.Modifier)
%     fprintf('%s released.\n',eventdata.Key{1})
    switch eventdata.Modifier{1}
        case {'shift','alt'}
            h.shift_on = 0;
            data.shift_on=h.shift_on;
            set(h.tg_DisplayMode,'UserData',data);
        case 'control'
            h.control_on = 0;
    end
end

guidata(hObject,h);


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonDownFcn(hObject, eventdata, h)
 CellSelectionButtonDownFcn(hObject,eventdata,h);


function SVMSelectionButtonDownFcn(hObject,eventdata,h)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pos=get(gca,'CurrentPoint');   
 
x = round(pos(1,1));
y  = round(pos(1,2));
x = min(max(1, round(x)), h.dat.res.Lx);
y = min(max(1, round(y)), h.dat.res.Ly);

h.dat.F.ichosen = h.dat.res.iclust1(y, x);
already_selected  = find(h.dat.F.ichosen_append==h.dat.F.ichosen);

if ~isempty(already_selected)
    h.dat.F.ichosen_append(already_selected)=[]; % unselect the same ROI if clicked again
else
    h.dat.F.ichosen_append = cat(1,h.dat.F.ichosen_append,h.dat.F.ichosen);
end


h=update_figure(h);
SelectionType =get(gcf,'SelectionType');

if h.control_on %% ----- with Ctrl- key ----- %%
    fprintf('Ctrl-')
    switch SelectionType
        case 'alt' % in fact, this is ctrl-left 
            fprintf('Left (add) \n');
%             h.dat.F.ichosen_append ;
        otherwise
            fprintf('%s\n',SelectionType);
    end
else 
    h.dat.F.ichosen_append =h.dat.F.ichosen; % reset

    switch SelectionType
        case 'alt' % Right click
            % flip currently selected unit
%             h.dat.cl.manual(h.dat.F.ichosen) = .5 - h.dat.cl.selected(h.dat.F.ichosen);
%             h = ApplyROIFilter(h);
%             h = buildLambdaValue(h);
        case 'open' % double click
            % unpin the manual selection on this cell
%             h.dat.cl.manual(h.dat.F.ichosen) = 0;
%             h = ApplyROIFilter(h);
%             h = buildLambdaValue(h);
        case 'extend' % shift left
%             h.dat.cl.redcell(h.dat.F.ichosen) = 1 -  h.dat.cl.redcell(h.dat.F.ichosen);
%             
%             if h.dat.cl.redcell(h.dat.F.ichosen)==1
%                 h.dat.cl.rands(h.dat.F.ichosen) = 0;
%             else
%                 h.dat.cl.rands(h.dat.F.ichosen) = h.dat.cl.rands_orig(h.dat.F.ichosen);
%             end
%             h.dat.graph.img1.H       = reshape(h.dat.cl.rands(h.dat.res.iclust), h.dat.res.Ly, h.dat.res.Lx);
%             h.dat.graph.img2.H       = reshape(h.dat.cl.rands(h.dat.res.iclust), h.dat.res.Ly, h.dat.res.Lx);
%             
%             if h.dat.cl.redcell(h.dat.F.ichosen)
%                 display('red')
%             else
%                 display('not red')
%             end
        case 'normal' % left click
%             fprintf('Left click\n')
            try delete(h.svm_editnow), catch, end
          
            xl = get(gca,'XLim');
            yl = get(gca,'YLim');
            pos_y = pos(1,2)
            pos_x = pos(1,1)
            npos_x = (pos_x-xl(1))/diff(xl);
            npos_y = (pos_y-yl(1))/diff(yl);
            dh = 0.03;
            dw = 0.03;
            pos=[npos_x-dw/2 1-npos_y-dh/2 dw dh];
%             pos=[0.5 0.5 0.03 0.03];
            CallbackFcn=@(hObject,eventdata)RoiGui_007('enter_svm_value', hObject,eventdata,h);            
            h.svm_editnow=uicontrol('Style', 'edit','Parent',h.uipanel_Left,'String', '0',...
                'Unit','normalized','Position',pos,...
                'UserData',h,...
                'Callback', CallbackFcn);
            waitfor(h.svm_editnow);
        otherwise
            SelectionType
    end
end


fprintf('\nclust=%d(iscell=%2.1f,manual=%2.1f), (x,y)=(%2d,%2d)\n',...
    h.dat.F.ichosen,...
    h.dat.cl.manual(h.dat.F.ichosen),...
    h.dat.cl.selected(h.dat.F.ichosen),x,y);
fprintf('Compactness=%2.1f, npix=%d\n', h.dat.cl.statTBL.Compactness(h.dat.F.ichosen),...
    h.dat.cl.statTBL.npix(h.dat.F.ichosen));
fprintf('Skew(F)=%3.3f, ',h.dat.cl.statTBL.skewF(h.dat.F.ichosen));
fprintf('Eccentricity=%3.3f, ',h.dat.cl.statTBL.Eccentricity(h.dat.F.ichosen));
fprintf('Solidity=%3.3f, ',h.dat.cl.statTBL.Solidity(h.dat.F.ichosen));
fprintf('meanV=%3.3f,',h.dat.cl.statTBL.V(h.dat.F.ichosen)/h.dat.cl.statTBL.npix(h.dat.F.ichosen));
fprintf('manual label=%3.3d\n',h.dat.cl.svm_class(h.dat.F.ichosen));

guidata(hObject,h);


function enter_svm_value(hObject,eventdata,h)
str=get(hObject,'String');
if ~ischar(str)
    delete(hObject);
    error('The cluster number must be a numeric!');
end
h=get(hObject,'UserData');

%<To Do> set manual id number of the cell.

delete(hObject);

function EllipseSelectionButtonDownFcn(hObject,eventdata,h)

pos=get(gca,'CurrentPoint');

x = round(pos(1,1));
y  = round(pos(1,2));
x = min(max(1, round(x)), h.dat.res.Lx);
y = min(max(1, round(y)), h.dat.res.Ly);


if isfield(h.dat,'ImEllipseH')
    live_ellipse = ishandle(h.dat.ImEllipseH);
    if all(live_ellipse==0)
        if isfield(h.dat,'ImEllipseCoveredArea')
            h.dat=rmfield(h.dat,'ImEllipseCoveredArea');
        end
    else
        h.dat.ImEllipseH=h.dat.ImEllipseH(live_ellipse);
        h.dat.ImEllipseCoveredArea = zeros(h.dat.res.Ly,h.dat.res.Lx);
        for ii=1:length(h.dat.ImEllipseH)
            h.dat.ImEllipseCoveredArea=h.dat.ImEllipseCoveredAre+h.dat.ImEllipseH(ii).createMask;
        end
    end
else
    
end

if isfield(h.dat,'ImEllipseCoveredArea')
    WithinImEllipse = h.dat.ImEllipseCoveredArea(y,x);
else
    WithinImEllipse=0;
end

if ~WithinImEllipse
    h=add_imellipse_Callback(hObject,eventdata,h);
end

h=pb_show_EllipticROI_Callback(hObject, eventdata, h);
guidata(hObject,h);

function CellSplitButtonDownFcn(hObject,eventdata,h)


pos=get(gca,'CurrentPoint');

x1 = round(pos(1,1));
y1  = round(pos(1,2));
x1 = min(max(1, round(x1)), h.dat.res.Lx);
y1 = min(max(1, round(y1)), h.dat.res.Ly);
hold on; plotH=plot(x1,y1,'ro');

% another center of k-means clustering

set(h.left_imageH,'HitTest','off'); % to avoid evoking the callback again 
x2 = ceil(ginput(1));
set(h.left_imageH,'HitTest','on');
y2 = x2(2);
x2 = x2(1);

delete(plotH);

if isfield(h.dat,'svd') && isfield(h.dat.svd,'U')
    % do nothing    
else
    ops=h.dat.ops;
  
    PlaneChString = sprintf('plane%d_ch%d',h.dat.ops.PlaneID,h.dat.ops.ChannelID);
    SVDFile = sprintf('%s/SVDroi_%s_%s_%s.mat', ops.ResultsSavePath, ...
        ops.mouse_name, ops.date, PlaneChString);
    if ~exist(SVDFile,'file')
        [LoadPath,~,~]=fileparts(h.dat.filename);
        SVDFile = sprintf('%s/SVDroi_%s_%s_%s.mat', LoadPath, ...
            ops.mouse_name, ops.date, PlaneChString);
    end
    fprintf('Loading svd file: %s\n', SVDFile);
    
    if exist(SVDFile,'file')    
        waitH = waitbar(0.2,sprintf('Loding %s',SVDFile),'Name','Loading SVD to calculate k-means..');
        h.dat.svd=   load(SVDFile, 'U','Sv');
        guidata(hObject,h);
        waitbar(1,waitH);
        delete(waitH);
    else
        error('SVD file %s NOT FOUND.',SVDFile);
    end
    
end

L=h.dat.res.iclust1;

target=L(y1,x1);
ind=find(L==target);
fprintf('target=%d\n ',target);

target2=L(y2,x2);
if target~=target2
    errordlg('different two region is selected');
end
% get the square region
[indY,indX]=ind2sub(size(L),ind);
center_ind=sub2ind(size(L),[y1,y2],[x1, x2]);

Val = L==target;

J=min(indY):max(indY);
I=min(indX):max(indX);

Sat=reshape(h.dat.res.probabilities,h.dat.res.Ly,h.dat.res.Lx);

% plot the principle components.
Utmp = reshape(h.dat.svd.U,[],length(h.dat.svd.Sv));
Utmp=Utmp.*sqrt(h.dat.svd.Sv'); % > 2016b
% U = bsxfun(@times, U, Sv'.^.5); < 20176b

Utmp=Utmp(ind,:);
Utmp2 = cat(2,zscore(indY),zscore(indX),zscore(Utmp,0,1)/5); % 
[coef,score] = pca(Utmp2);

% Col = {'r','g','b'}; 
% grid on;
Nlabel = 2;

ind1=find(ind==center_ind(1));
ind2=find(ind==center_ind(2));
idx = kmeans(Utmp2,Nlabel,'Start',Utmp2([ind1,ind2],:));

% make the divided image 

Col = [0.9,0.4,0.6];
Hue = h.dat.graph.img1.H;
for ii=1:Nlabel
    original_ind = ind(idx==ii);
    Hue(original_ind)=Col(ii);
end

IMG=hsv2rgb(cat(3,Hue,Sat,Val));

SplitOrNot = ROISplit_GUI('title','ROI split',...
    'CData',IMG(J,I,:),...
    'PlotData',score,...
    'PlotData_Class',idx)

unique_roi = unique(h.dat.res.iclust1);
nROI = max(unique_roi);
switch SplitOrNot
    case 'Split'
        
        for ii=1:Nlabel
            BW=zeros(size(h.dat.res.iclust1));
            BW(ind(idx==ii))=1;
            CC = bwconncomp(BW);
            [~,max_id] = max(cellfun(@length,CC.PixelIdxList));
            % here, pickup the largest ROI as sub-divided area, and other 
%             tiny fragments are left with original iclust number. 
            nROI= nROI+1;
            h.dat.res.iclust1(CC.PixelIdxList{max_id})=nROI;
            tmp= get_stat_from_iclust(h.dat.res.iclust1==nROI,h.dat.res.M);
            tmp.skewF=inf;
            tmp.SNratio=1;
            h.dat.stat(nROI)=tmp;
        end
        %exclude original ROI
     
        h.dat.cl.manual(target)=-0.5;
        h.dat.cl.selected(target)=0; % set the original unselected. 
        
        h=update_cl(h);
        fprintf('iclust=%d is split into %d and %d\n',target,nROI-1,nROI);
        % <To Do>
        % return to normal mode
        
    case 'No, I change my mind'
        
    otherwise
        error('Unknown answer %s',SplitOrNot)
end

% come back to display mode        

% set(h.left_imageH,'Hittest','on');
% h.ModeSelectionButtonGroup.SelectedObject=h.tg_DisplayMode;

h=redraw_figure(h);
h=update_figure(h);
guidata(hObject,h);

 % cell selection function
function CellSelectionButtonDownFcn(hObject,eventdata,h)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pos=get(gca,'CurrentPoint');

x = round(pos(1,1));
y  = round(pos(1,2));
x = min(max(1, round(x)), h.dat.res.Lx);
y = min(max(1, round(y)), h.dat.res.Ly);

if isempty(h.dat.F.ichosen_append), h.dat.F.ichosen_append=[];end

if  h.dat.graph.img0.BackGround(y,x)
    disp('Background clicked');
    return;
end

h.dat.F.ichosen = h.dat.res.iclust1(y, x);
already_selected  = find(h.dat.F.ichosen_append==h.dat.F.ichosen);

if ~isempty(already_selected)
    h.dat.F.ichosen_append(already_selected)=[]; % unselect the same ROI if clicked again
else
    h.dat.F.ichosen_append = cat(1,h.dat.F.ichosen_append,h.dat.F.ichosen);
end

SelectionType =get(gcf,'SelectionType');

if h.control_on %% ----- with Ctrl- key ----- %%
    fprintf('Ctrl-')
    switch SelectionType
        case 'alt' % in fact, this is ctrl-left 
            fprintf('Left (add) \n');
%             h.dat.F.ichosen_append ;
        otherwise
            fprintf('%s\n',SelectionType);
    end
else 
    h.dat.F.ichosen_append =h.dat.F.ichosen; % reset

    switch SelectionType
        case 'alt' % Right click
            % flip currently selected unit
            h.dat.cl.manual(h.dat.F.ichosen) = .5 - h.dat.cl.selected(h.dat.F.ichosen);
            h = ApplyROIFilter(h);
            h = buildLambdaValue(h);
        case 'open' % double click
            % unpin the manual selection on this cell
            h.dat.cl.manual(h.dat.F.ichosen) = 0;
            h = ApplyROIFilter(h);
            h = buildLambdaValue(h);
        case 'extend' % shift left
            h.dat.cl.redcell(h.dat.F.ichosen) = 1 -  h.dat.cl.redcell(h.dat.F.ichosen);
            
            if h.dat.cl.redcell(h.dat.F.ichosen)==1
                h.dat.cl.rands(h.dat.F.ichosen) = 0;
            else
                h.dat.cl.rands(h.dat.F.ichosen) = h.dat.cl.rands_orig(h.dat.F.ichosen);
            end
            h.dat.graph.img1.H       = reshape(h.dat.cl.rands(h.dat.res.iclust1), h.dat.res.Ly, h.dat.res.Lx);
            h.dat.graph.img2.H       = reshape(h.dat.cl.rands(h.dat.res.iclust1), h.dat.res.Ly, h.dat.res.Lx);
            
            if h.dat.cl.redcell(h.dat.F.ichosen)
                display('red')
            else
                display('not red')
            end
        case 'normal' % left click
%             fprintf('Left click\n')
        otherwise
            SelectionType
    end
end



fprintf('\nclust=%d(iscell=%2.1f,manual=%2.1f), (x,y)=(%2d,%2d)\n',...
    h.dat.F.ichosen,...
    h.dat.cl.selected(h.dat.F.ichosen),...
    h.dat.cl.manual(h.dat.F.ichosen),...
    x,y);
try
    fprintf('Compactness=%2.1f, npix=%d\n', h.dat.cl.statTBL.Compactness(h.dat.F.ichosen),...
        h.dat.cl.statTBL.npix(h.dat.F.ichosen));
catch
    fprintf('Compactness=%2.1f, npix=%d\n', h.dat.cl.statTBL.Compactness{h.dat.F.ichosen},...
        h.dat.cl.statTBL.npix{h.dat.F.ichosen});
end
fprintf('Skew(F)=%3.3f, ',h.dat.cl.statTBL.skewF(h.dat.F.ichosen));
try fprintf('Eccentricity=%3.3f, ',h.dat.cl.statTBL.Eccentricity(h.dat.F.ichosen));
    fprintf('Solidity=%3.3f, ',h.dat.cl.statTBL.Solidity(h.dat.F.ichosen));
    fprintf('meanV=%3.3f,',h.dat.cl.statTBL.V(h.dat.F.ichosen)/h.dat.cl.statTBL.npix(h.dat.F.ichosen));
catch
    fprintf('Eccentricity=%3.3f, ',h.dat.cl.statTBL.Eccentricity{h.dat.F.ichosen}); 
    fprintf('Solidity=%3.3f, ',h.dat.cl.statTBL.Solidity{h.dat.F.ichosen});
    fprintf('meanV=%3.3f,',h.dat.cl.statTBL.V{h.dat.F.ichosen}/h.dat.cl.statTBL.npix{h.dat.F.ichosen});
end

fprintf('manual label=%d\n',h.dat.cl.svm_class(h.dat.F.ichosen));


h=update_figure(h);
guidata(hObject,h);

function h=update_figure(h)
Mode = h.ModeSelectionButtonGroup.SelectedObject.Tag;

switch Mode
    case {'tg_DisplayMode','tg_ROISplitMode','tg_ROIMergeMode'}
        PlotHue = h.dat.cl.rands(h.dat.F.ichosen_append);
        PlotColor =permute(cat(2,PlotHue(:),ones(length(PlotHue),2)),[1,3,2]);
        PlotColor = hsv2rgb(PlotColor);
        
        h=redraw_fluorescence_multi(h,PlotColor);
        
        h = buildSat(h);
        h = buildLambdaValue(h);
        h = buildHue(h);
        
        %fprintf('calculating cell mask\n')
%         h.dat.cellmask=zeros(h.dat.res.Ly, h.dat.res.Lx);
%         for ii=1:length(h.dat.F.ichosen_append)
% %             h.dat.cellmask(h.dat.res.iclust==h.dat.F.ichosen_append(ii))=1;
%             h.dat.cellmask(h.dat.stat(h.dat.F.ichosen_append(ii)).ipix_edge)=1;
%         end
%         h.dat.cellmask = edge(h.dat.cellmask);
        h=redraw_figure(h);
        
    case 'tg_EllipseMode'
%         h=add_imellipse_Callback(hObject, eventdata, h);
%     case 'pb_ROIMergeMode'
%         h=redraw_corrimg(h);
    case 'tb_SVMMode'
        PlotHue = h.dat.cl.rands(h.dat.F.ichosen_append);
        PlotColor =permute(cat(2,PlotHue(:),ones(length(PlotHue),2)),[1,3,2]);
        PlotColor = hsv2rgb(PlotColor);
        h=redraw_fluorescence_multi(h,PlotColor);
        
        h = buildHue(h);
        h = buildSat(h);
        h = buildLambdaValue(h);
    
        
        %fprintf('calculating cell mask\n')
%         h.dat.cellmask=zeros(h.dat.res.Ly, h.dat.res.Lx);
%         for ii=1:length(h.dat.F.ichosen_append)
%             h.dat.cellmask(h.dat.res.iclust==h.dat.F.ichosen_append(ii))=1;
%         end
%         h.dat.cellmask = edge(h.dat.cellmask);
        h=redraw_figure(h);
    case 'tb_LabelMode'
        
        ROIColor=h.popup_ROI_ColorSelection.String{h.popup_ROI_ColorSelection.Value};
        switch ROIColor
            case 'Green'
                PlotHue =0.333;
            case 'Red'
                PlotHue=0;
            case 'Blue'
                PlotHue=0.667;
            otherwise
                error('Unknown color %s',ROIColor);
        end
        PlotHue = PlotHue*ones(length(h.dat.F.ichosen_append),1);
        
        PlotColor =permute(cat(2,PlotHue(:),ones(length(PlotHue),2)),[1,3,2]);
        PlotColor = hsv2rgb(PlotColor);
        
        h=redraw_fluorescence_multi(h,PlotColor);
        
        h = buildSat(h);
        h = buildLambdaValue(h);
        h = buildHue(h);
        
        h=redraw_figure(h);
        
           
    otherwise
        error('Unknown Tag %s',Mode);
end



function ROISelectionButtonDownFcn(hObject,eventdata,h);

pos=get(gca,'CurrentPoint')

% eventdata
% z = round(eventdata.Source.CurrentAxes.CurrentPoint(1,:));
x = round(pos(1,1));
y  = round(pos(1,2));

% disp(eventdata.Source.SelectionType)
% keyboard;
% manual selection of units
x = min(max(1, round(x)), h.dat.res.Lx);
y = min(max(1, round(y)), h.dat.res.Ly);

SelectionType =get(gcf,'SelectionType')

switch SelectionType
    case 'normal'
        fprintf('Frame=%d\n',x)
    case 'alt' % Right click
     
    case 'open' % double click
      
    case 'extend' % shift left
     
end




% --- Executes on button press in pb_proccdisplay.
function h=pb_proccdisplay_Callback(hObject, eventdata, h)

h.dat.procmap = 1 -  h.dat.procmap;
if (h.dat.procmap)
    disp('ProcMap On (Normalized image)');
else
    disp('ProcMap Off (Normalized Image)');
end

if h.dat.graph.map>1
    h=redraw_meanimg(h);
end

% if h.dat.graph.map < h.dat.graph.maxmap
%     h.dat.graph.map = h.dat.graph.map + 1;
%     h=redraw_meanimg(h);
% else
%     h.dat.graph.map = 1;
%     h=redraw_figure(h);
% end
guidata(hObject,h);


% --- Executes on button press in pb_rois_display.
function h=pb_rois_display_Callback(hObject, eventdata, h)
h.dat.graph.map = 1;
h=redraw_figure(h);
guidata(hObject,h);


% --- Executes on button press in pb_meandisplay.
function h=pb_meandisplay_Callback(hObject, eventdata, h)
h.dat.graph.map = 2;
h=redraw_meanimg(h);
guidata(hObject,h);

% --- Executes on button press in pb_Rede_display.
function pb_Rede_display_Callback(hObject, eventdata, h)
 h.dat.graph.map = 3;
h=redraw_meanimg(h);
guidata(hObject,h);

% RED CORRECTED BUTTON
function pb_redcorr_display_Callback(hObject, eventdata, h)
h.dat.graph.map = 4;
h=redraw_meanimg(h);
guidata(hObject,h);






% --- Executes on button press in rb_showneuropil.
function rb_showneuropil_Callback(hObject, eventdata, h)
% h.dat.plot_neu = 1 - h.dat.plot_neu;
h=redraw_fluorescence(h);
guidata(hObject,h);



% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over pb_proccdisplay.
function pb_proccdisplay_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to pb_proccdisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function h=redraw_fluorescence_multi(h,varargin)

if nargin>=2
    PlotColor = varargin{1};
else
    PlotColor = repmat([0,0,1],length(h.plot_fluorescence_F),1); 
end

[NN NT] = size(h.dat.F.trace);
h.dat.F.ichosen_append(h.dat.F.ichosen_append>NN) =[];
if isempty(h.dat.F.ichosen_append)
    if isempty(h.dat.F.ichosen_append)
        x_trace = [1 NT];
        y_trace=  [0 0];
    else
        x_trace = [1:NT];
         if get(h.rb_showneuropil,'Value')
             % neuropil subtracted 
              y_trace=my_conv_local(double(h.dat.F.truetrace(h.dat.F.ichosen,:)), 3);
         else
             % raw fluorescence. 
              y_trace=my_conv_local(double(h.dat.F.trace(h.dat.F.ichosen,:)), 3);
         end
    end
else
    x_trace = [1:NT]; %repmat([1:NT NaN],size(h.dat.F.ichosen_append,1),1)';
    if get(h.rb_showneuropil,'Value')
      y_trace=zscore(my_conv_local(double(h.dat.F.truetrace(h.dat.F.ichosen_append,:)), 3),0,2);
    else
       y_trace=zscore(my_conv_local(double(h.dat.F.trace(h.dat.F.ichosen_append,:)), 3),0,2);
    end
    
 
    y_trace = bsxfun(@plus, y_trace,2*[0:length(h.dat.F.ichosen_append)-1]');
%     y_trace = cat(2,y_trace,NaN(size(y_trace,1),1))';
%     
%     x_trace=x_trace(:);
%     y_trace=y_trace(:);
end

for ii=1:length(h.plot_fluorescence_F)
    if ii<=length(h.dat.F.ichosen_append)
        set(h.plot_fluorescence_F(ii),'XData',x_trace,'YData',y_trace(ii,:),'Color',PlotColor(ii,:));
    else
        set(h.plot_fluorescence_F(ii),'XData',NaN,'YData',NaN,'Color',[0 0 1]);
    end
end
 max_y = nanmax(y_trace(:));
 min_y = nanmin(y_trace(:));
       
%      if get(h.rb_showneuropil,'Value') && isfield(h.dat.F, 'neurop')
%          y_neurop=my_conv_local(double(h.dat.F.neurop(h.dat.F.ichosen,:)), 3);
%          x = [ones(length(y_neurop),1), y_neurop'];
%          %      b=regress(y_trace,x);
%          b = x'*x\x'*y_trace';
%          Coef = 0.5;
%          y_neurop = Coef*b(2)*y_neurop;
%          
%          set(h.plot_fluorescence_baseF,'XData',1:NT,'YData',y_neurop-mean(y_neurop)+mean(y_trace)/2,'Color','c','LineStyle','--');
%          max_y = max([max_y; y_neurop(:)]);
%          min_y = min([min_y; y_neurop(:)]);
%       
%          newydata=y_trace-Coef*b(2)*y_neurop+mean(y_trace)*1.1;
%          set(h.plot_fluorescence_FCell,'XData',1:NT,'YData',newydata,'Color','r');
%      else
%          set(h.plot_fluorescence_FCell,'XData',NaN,'YData',NaN,'Color','r');
%          set(h.plot_fluorescence_baseF,'XData',NaN,'YData',NaN);
%      end

 dy = max_y-min_y;
 if dy==0,    dy=1;end
 yl(1) = min_y-0.1*dy;
 yl(2)=  max_y+0.1*dy;
 
 set(h.axes_fluorescence,'XLim',[0 NT],'YLim',yl);
 set(h.plot_fluorescence_baseF,'XData',NaN,'YData',NaN);
  set(h.plot_fluorescence_FCell,'XData',NaN,'YData',NaN,'Color','r');
guidata(h.axes_fluorescence,h);

function h=redraw_fluorescence(h)
% hold off
[NN NT] = size(h.dat.F.trace);
% plot(my_conv_local(medfilt1(double(h.dat.F.trace(h.dat.F.ichosen,:)), 3), 3));
% ydata=my_conv_local(medfilt1(double(h.dat.F.trace(h.dat.F.ichosen,:)), 3), 3);
% ydata=double(h.dat.F.trace(h.dat.F.ichosen,:));
 y_trace=my_conv_local(double(h.dat.F.trace(h.dat.F.ichosen,:)), 3);
 
%  try 
%      delete(get(h.axes_fluorescence,'Children'))
%  end
set(h.plot_fluorescence_F,'XData',1:NT,'YData',y_trace,'Color','b');

% set(h.plot_fluorescence_F,'XData',NaN,'YData',NaN,'Color','c');
max_y = max(y_trace);
min_y = min(y_trace);
% dy = max_y-min_y;
% if dy==0,    dy=1;end
% yl(1) = min_y-0.1*dy;
% yl(2)=  max_y+0.1*dy;
% 
% 
% set(h.axes_fluorescence,'XLim',[0 NT],'YLim',yl);


     
% ydata
% axis tight;

% fprintf('Var(F)=%3.3f\n',var(ydata));
fprintf('Skew(F)=%3.3f\n',h.dat.cl.statTBL.skewF(h.dat.F.ichosen));

if get(h.rb_showneuropil,'Value')
    if isfield(h.dat.F, 'neurop')
      
        y_neurop=my_conv_local(double(h.dat.F.neurop(h.dat.F.ichosen,:)), 3);

%         x = [ones(length(y_neurop),1), y_neurop'];
%         % b=regress(y_trace,x);
%         b = x'*x\x'*y_trace';
%         Coef = 0.6;
%         newydata=y_trace-Coef*b(2)*y_neurop-b(1);
%         newydata=my_conv_local(newydata',3)+mean(y_trace);


        set(h.plot_fluorescence_baseF,'XData',1:NT,'YData',y_neurop,'Color','c','LineStyle','--');
%         max_y = max([max_y; y_neurop(:)]);
%         min_y = min([min_y; y_neurop(:)]);
    else
        set(h.plot_fluorescence_baseF,'XData',NaN,'YData',NaN);

    end
else
    y_neurop = [];
    set(h.plot_fluorescence_baseF,'XData',NaN,'YData',NaN);
end

if isempty(y_neurop)
    newydata = y_trace;
else
    Coef='1/Regress';
    newydata = subtract_neurop(y_trace,y_neurop,Coef);
    newydata=my_conv_local(newydata',3)+mean(y_trace);
end

set(h.plot_fluorescence_FCell,'XData',1:NT,'YData',newydata,'Color','r');
max_y = max([max_y;newydata(:)]);
min_y = min([min_y;newydata(:)]);
dy = max_y-min_y;
if dy==0,    dy=1;end
yl(1) = min_y-0.1*dy;
yl(2)=  max_y+0.1*dy;

set(h.axes_fluorescence,'XLim',[0 NT],'YLim',yl);

guidata(h.axes_fluorescence,h);
% plot([0 NT], [0 0], 'k', 'Linewidth', 2)
% axis off

function h=redraw_figure(h)

if h.dat.display_select % show SELECTED one  
    I = hsv2rgb(cat(3, h.dat.graph.img1.H, h.dat.graph.img1.Sat, h.dat.graph.img1.V));
else % show UN-selected one
    I = hsv2rgb(cat(3, h.dat.graph.img2.H, h.dat.graph.img2.Sat, h.dat.graph.img2.V));
end

I = min(I, 1);
% axes(h.axes_left); imagesc(I);
set(h.left_imageH,'CData',I);
set(h.axes_left,'XLim',[h.dat.QV.xlim],'YLim',[h.dat.QV.ylim]);



    function h=redraw_meanimg(h)

if h.dat.procmap
    I = h.dat.graph.mimg_proc(:,:,h.dat.graph.map);
else    
    I = h.dat.graph.mimg(:,:,h.dat.graph.map);
end

mu = median(I(:));
sd1 = mean(abs(I(I<mu) - mu));
sd2 = mean(abs(I(I>mu) - mu));

% axes(h.axes_left); imagesc(I, mu + 5*[-sd1 sd2]);
set(h.left_imageH,'CData',I);
% xlim([h.dat.QV.xlim]); ylim([h.dat.QV.ylim]);
set(h.axes_left,'XLim',[h.dat.QV.xlim],'YLim',[h.dat.QV.ylim],'CLim',mu + [-2.5*sd1 4*sd2]);

% axes(h.axes_right); imagesc(I, mu + 5*[-sd1 sd2]);
% xlim([h.dat.QV.xlim]); ylim([h.dat.QV.ylim]);
% set(h.left_imageH,'CData',I);
% xlim([h.dat.QV.xlim]); ylim([h.dat.QV.ylim]);
% set(h.axes_right,'XLim',[h.dat.QV.xlim],'YLim',[h.dat.QV.ylim],'CLim',mu + [-2.5*sd1 4*sd2]);


drawnow

function h=show_corr_roi(h)
%%
IDX=h.dat.cl.IDX;
unique_cluster = unique(IDX);
unique_cluster(unique_cluster==0)=[];



% roi_samecell=find(h.dat.cl.C_of_zF(h.dat.F.ichosen,:)>threshold);


Vmap = h.dat.graph.img1.V; % value, 1 is brighter


figh=myfigure('CorrROI');clf;% set(figh,'Position',);
for ii=unique_cluster(:)'
    subplot(2,1,1);
    roi_samecell = find(IDX==ii);
    Hmap = zeros(1,length(h.dat.res.iclust1)); % hue color map, make the same cell same color.
    Hmap(roi_samecell)=1;
    Hmap       = reshape(Hmap(h.dat.res.iclust1), h.dat.res.Ly, h.dat.res.Lx);
    Smap = zeros(size(h.dat.graph.img1.Sat)); % saturation, 1 is more pure color, 0 is white.
    for jj=roi_samecell(:)'
        Smap(h.dat.res.iclust1==jj)=1;
    end

    I = hsv2rgb(cat(3, Hmap, Smap, Vmap));
    I = min(I,1);
    imagesc(I);
    title_txt = sprintf('Cluster=%d ',ii);
    title(title_txt)
    
    subplot(2,1,2);
    F = zscore(h.dat.F.trace(roi_samecell,:),0,2);
    F=my_conv_local(F, 3);
    plot(1:size(F,2),F);
    
    pause
end

function h=redraw_corrimg(h,varargin)

% dbscan to find correlated cluster.
% epsilon = str2double(  get(h.edit_CorrROI_epsilon ,'String'));
% MinPts = str2double(  get(h.edit_CorrROI_MinPts,'String'));
% IDX=dbscan_D(1-h.dat.cl.C_of_zF,epsilon,MinPts);
IDX = h.dat.cl.IDX; %recalc_IDX(h);

if IDX(h.dat.F.ichosen)~=0
    roi_samecell = find(IDX==IDX(h.dat.F.ichosen));
else
    roi_samecell = h.dat.F.ichosen;
end

selected_region = edge(h.dat.res.iclust1==h.dat.F.ichosen);

Smap = zeros(size(h.dat.graph.img1.Sat)); % saturation, 1 is more pure color, 0 is white.
% color the selected cell, and let the others remain white. 
for ii=1:length(roi_samecell)
    Smap(h.dat.res.iclust1==roi_samecell(ii))=1;
end

Hmap      = reshape(h.dat.cl.rands(h.dat.res.iclust1), h.dat.res.Ly, h.dat.res.Lx);
% 
% % Hmap = zeros(1,length(h.dat.res.iclust1)); % hue color map, make the same cell same color.
% % Hmap(roi_samecell)=0.6; % blue
% % Hmap=h.dat.cl.rands;
% % Hmap       = reshape(Hmap(h.dat.res.iclust1), h.dat.res.Ly, h.dat.res.Lx);


Vmap = h.dat.res.M0; 
%h.dat.graph.img1.V+selected_region; % value, 1 is brighter



I = hsv2rgb(cat(3, Hmap, Smap, Vmap));
I = min(I, 1);
% axes(h.axes_left); imagesc(I);
set(h.left_imageH,'CData',I);
set(h.axes_left,'XLim',[h.dat.QV.xlim],'YLim',[h.dat.QV.ylim]);
h.dat.F.ichosen_append = roi_samecell;


PlotHue = h.dat.cl.rands(h.dat.F.ichosen_append);
PlotColor =permute(cat(2,PlotHue(:),ones(length(PlotHue),2)),[1,3,2]);
PlotColor = hsv2rgb(PlotColor);

h=redraw_fluorescence_multi(h,PlotColor);




function h=redraw_movie(hObject,eventdata,h)

if ~isfield(h.dat,'reg_data')
    ButtonName = questdlg('Regx5 file not loaded yet. Load? ', ...
        'Regx5 load question', ...
        'YES', 'NO', 'YES');
    switch ButtonName
        case 'YES'
            h=local_loadmovie(hObject, eventdata, h);
        case 'NO'
            disp('Use cancelled');
            return;
    end
end
sd1=0;
sd2=1025;
% RightH = h.left_imageH;
% set(RightH,'CData',h.dat.reg_data.get_single_slice(1));

set(h.left_imageH,'CData',h.dat.reg_data.get_single_slice(1));
set(h.axes_left,'CLim',[sd1 sd2],'UserData',1);


function previous_frame(hObject,eventdata,h)

%// necessary to check if the timer is already running
%// otherwise the automatic key repetition tries to start
%// the timer multiple time, which produces an error
if strcmp(h.movie_timer.Running,'off')
    userdata=get(h.axes_left,'UserData');
    userdata(2)=-1; % role back.
    set(h.axes_left,'UserData',userdata);
    start(h.movie_timer)
end

function goto_frame(hObject,eventdata,h,index)
% to translate movie index to fast reg-movie index.
index = round(index/h.dat.reg_decimation);

userdata=get(h.axes_left,'UserData');
userdata(1)=index; % role back.
set(h.axes_left,'UserData',userdata);

set(h.left_imageH,'CData',h.dat.reg_data.get_single_slice(index));

set(h.uipanel_bottom,'Title',sprintf('(%d/%d)',index,h.dat.reg_data.TotalNSeries));
set(h.plot_fluorescence_timing_bar,'XData',index*h.dat.reg_decimation*[1 1],...
    'YData',get(h.axes_fluorescence,'YLim'));
% 
% 
drawnow;


function next_frame(hObject,eventdata,h)

if strcmp(h.movie_timer.Running,'off')
    userdata=get(h.axes_left,'UserData');
    userdata(2)=1; % role forward.
    set(h.axes_left,'UserData',userdata);
    start(h.movie_timer)
end


function mimg = normalize_image(mimg)

dF = mimg - my_conv2(mimg, 2, [1 2]);
mimg = dF ./ my_conv2(abs(dF), 4, [1 2]);

% mimg = mimg - my_conv2(mimg, 5, [1 2]);
% mimg = mimg ./ my_conv2(abs(mimg), 10, [1 2]);

function h = get_parent_stats(h)

h.dat.stat(1).V = [];
for i = 1:length(h.dat.stat)
    if ~isfield(h.dat.stat(i), 'V') || isempty(h.dat.stat(i).V)
        h.dat.stat(i).region;
        h.dat.stat(i).V         = sum([h.dat.stat(i).region.V]);
        h.dat.stat(i).mrs = min(h.dat.stat(i).mrs, 1e4);
       % this is a parent region 
        h.dat.stat(i).parent    = i;
        h.dat.stat(i).VperPix   = h.dat.stat(i).V/h.dat.stat(i).npix;
        h.dat.stat(i).npix_res  = numel(h.dat.stat(i).ipix);
        h.dat.stat(i).nregions  = numel(h.dat.stat(i).region);
    end
end

nreg = [h.dat.stat.nregions];
npix_res = [h.dat.stat.npix_res];
npix = [h.dat.stat.npix];

VperPix = [h.dat.stat.VperPix];
mrs = [h.dat.stat.mrs]./[h.dat.stat.mrs0];
iparent = [h.dat.stat.parent];

h.dat.cl.nreg       = nreg(iparent);
h.dat.cl.npix_res   = npix_res(iparent);
h.dat.cl.npix_par   = npix(iparent);
h.dat.cl.VperPix    = VperPix(iparent);
h.dat.cl.mrs_parent = mrs(iparent);

function h = ApplyROIFilter(h)

contents=get(h.popup_FilterSelection,'String');
SelectedFilter = contents{get(h.popup_FilterSelection,'Value')};


switch SelectedFilter
    case 'Param filter'
        h.dat.cl.selected = h.dat.cl.statTBL.Compactness <h.dat.res.Mrs_thresh & ...
            h.dat.cl.statTBL.npix <h.dat.cl.npix_high & ...
            h.dat.cl.statTBL.npix >h.dat.cl.npix_low;
        
        % h.dat.cl.selected = h.dat.cl.selected  & (h.dat.cl.nreg      <h.dat.cl.nreg_max);
        % h.dat.cl.selected = h.dat.cl.selected  & (h.dat.cl.npix_par  <h.dat.cl.npix_par_max);
        % h.dat.cl.selected = h.dat.cl.selected  & (h.dat.cl.npix_res  <h.dat.cl.npix_res_max);
        % h.dat.cl.selected = h.dat.cl.selected  & (h.dat.cl.mrs_parent<h.dat.cl.mrs_parent_max);
        % h.dat.cl.selected = h.dat.cl.selected  & (h.dat.cl.VperPix   >h.dat.cl.VperPix_min);
        h.dat.cl.skewF_low = str2double(get(h.edit_minSkewF,'String'));
        h.dat.cl.EccentricityMax = str2double(get(h.edit_EccentricityMax,'String'));
        h.dat.cl.SolidityMin = str2double(get(h.edit_SolidityMin,'String'));
        
        h.dat.cl.selected = h.dat.cl.selected  & (h.dat.cl.statTBL.skewF     >h.dat.cl.skewF_low);
        h.dat.cl.selected = h.dat.cl.selected  & (h.dat.cl.statTBL.Eccentricity     <h.dat.cl.EccentricityMax);
        h.dat.cl.selected = h.dat.cl.selected  & (h.dat.cl.statTBL.Solidity     >h.dat.cl.SolidityMin);
        
       
    case 'SVM filter'
        h.dat.cl.selected = zeros(size(h.dat.cl.statTBL,1),1);
        
        if (h.cb_SVM_CellBody.Value)
            h.dat.cl.selected(h.dat.cl.svm_predicted_type==1)=1;
        end
        
        if (h.cb_SVM_Dendrite.Value)
            h.dat.cl.selected(h.dat.cl.svm_predicted_type==2)=1;
        end
        
        if (h.cb_SVM_Axon.Value)
            h.dat.cl.selected(h.dat.cl.svm_predicted_type==3)=1;
        end
        
         h.dat.cl.selected(h.dat.cl.svm_predicted_type==9)=0;
        
%          % remove the background.  -> move to init 
         backgroundID= unique(h.dat.res.iclust(h.dat.res.probabilities==0));
         h.dat.cl.manual(backgroundID)=-0.5;
         
    case 'none'
        
    otherwise
        error('Unknown filter %s',SelectedFilter);
end

% h.dat.cl.selected(max(h.dat.res.iclust1(:)))=0; % not anymore
h.dat.cl.selected = double(h.dat.cl.selected);

% always manual selections can overwrite. 
h.dat.cl.selected(h.dat.cl.manual>1e-3) = 1;
h.dat.cl.selected(h.dat.cl.manual<-1e-3) = 0;
h.dat.cl.k1 = reshape(h.dat.cl.selected(h.dat.res.iclust1), h.dat.res.Ly, h.dat.res.Lx);

function h = setOriginalThresh(h)
h.dat.res.Mrs_thresh = h.dat.res.Mrs_thresh_orig;
h.dat.cl.npix_low = h.dat.cl.npix_low_orig;
h.dat.cl.npix_high = h.dat.cl.npix_high_orig;

function h = buildHue(h)
Tag=h.ModeSelectionButtonGroup.SelectedObject.Tag;
h.dat.graph.img1.H       = reshape(h.dat.cl.rands(h.dat.res.iclust1), h.dat.res.Ly, h.dat.res.Lx);
h.dat.graph.img2.H       = h.dat.graph.img1.H;
        
switch Tag
    case {'tg_DisplayMode','tg_EllipseMode','tg_ROIMergeMode','tg_ROISplitMode'}    
     % do nothing
    case 'tb_SVMMode'
        
       for ii=1:h.dat.ops.Nk
            %             h.dat.cellmask(h.dat.res.iclust1==h.dat.F.ichosen_append(ii))=1;
            edge_pix = h.dat.stat(ii).ipix_edge;
            type_index = find([h.dat.cl.type_number_table{:,1}]==h.dat.cl.svm_class(ii));
            HSV=rgb2hsv(h.dat.cl.type_number_table{type_index,3}); 
            if type_index~=1
                h.dat.graph.img1.H(edge_pix)=HSV(1);
                h.dat.graph.img2.H(edge_pix)=HSV(1);
            end
       end
    case 'tb_LabelMode'
        % overwrite the hue value 
        
         ROIColor=h.popup_ROI_ColorSelection.String{h.popup_ROI_ColorSelection.Value};
        switch ROIColor
            case 'Green'
                PlotHue =0.333;
            case 'Red'
                PlotHue=0;
            case 'Blue'
                PlotHue=0.667;
            otherwise
                error('Unknown color %s',ROIColor);
        end
        PlotHue=PlotHue*ones(size(h.dat.cl.rands));
        
        h.dat.graph.img1.H       = reshape(PlotHue(h.dat.res.iclust1), h.dat.res.Ly, h.dat.res.Lx);
        h.dat.graph.img2.H       = h.dat.graph.img1.H;

    otherwise
        error('Unknown Tag %s',Tag);
end


function h = buildSat(h) % control gray to color mode. 0 is gray. 1 is color.

Tag=h.ModeSelectionButtonGroup.SelectedObject.Tag;

if isfield(h.dat.res,'probabilities')
    h.dat.graph.img0.BackGround = reshape(h.dat.res.probabilities==0,h.dat.res.Ly,h.dat.res.Lx);
else
    h.dat.graph.img0.BackGround = h.dat.res.iclust1==max(h.dat.res.iclust1(:));
end


% Sat=~h.dat.graph.img0.BackGround; % make the background gray color.
 Sat = reshape(sqrt(h.dat.res.probabilities),h.dat.res.Ly,h.dat.res.Lx);
switch Tag
    case {'tg_DisplayMode','tg_EllipseMode','tg_ROIMergeMode','tg_ROISplitMode','tb_LabelMode'}
              
        % make selected cell white fringed.
        h.dat.cellmask=zeros(h.dat.res.Ly, h.dat.res.Lx);
        for ii=1:length(h.dat.F.ichosen_append)
%             h.dat.cellmask(h.dat.res.iclust==h.dat.F.ichosen_append(ii))=1;
            edge_pix = h.dat.stat(h.dat.F.ichosen_append(ii)).ipix_edge;
            h.dat.cellmask(edge_pix)=1;
            Sat(edge_pix)=0;
        end
           
        h.dat.graph.img1.Sat     = Sat;
        h.dat.graph.img2.Sat     = Sat;
   
    case 'tb_SVMMode'
      
        for ii=1:h.dat.ops.Nk
            %             h.dat.cellmask(h.dat.res.iclust==h.dat.F.ichosen_append(ii))=1;
            edge_pix = h.dat.stat(ii).ipix_edge;
            type_index = find([h.dat.cl.type_number_table{:,1}]==h.dat.cl.svm_class(ii));
            HSV=rgb2hsv(h.dat.cl.type_number_table{type_index,3}); 
            if type_index~=1
                Sat(edge_pix)=HSV(2);
            end
        end
        
        h.dat.cellmask=zeros(h.dat.res.Ly, h.dat.res.Lx);
        for ii=1:length(h.dat.F.ichosen_append)
            %             h.dat.cellmask(h.dat.res.iclust==h.dat.F.ichosen_append(ii))=1;
            edge_pix = h.dat.stat(h.dat.F.ichosen_append(ii)).ipix_edge;
            h.dat.cellmask(edge_pix)=1;
            Sat(edge_pix)=0;
        end
        
        h.dat.graph.img1.Sat     = Sat;
        h.dat.graph.img2.Sat     = Sat;
        
    otherwise
        error('Unkown mode %s',Tag);
end


function h = buildLambdaValue(h) % V of HSV color; it controls Brightness. 
% h.dat.graph.img1.V       = h.dat.graph.img0.V .* h.dat.cl.k1 .* (1-h.dat.cl.excluded_pixels)...
%     .*(1-h.dat.cl.excluded_regions).*(1-h.dat.cl.excl_pix_perc).*h.dat.cl.topregion;
% h.dat.graph.img2.V       = h.dat.graph.img0.V .* (~h.dat.cl.k1);
Tag=h.ModeSelectionButtonGroup.SelectedObject.Tag;

BGval = get(h.slider_BGLevel,'Value');
BG=BGval*h.dat.graph.img0.BackGround.*h.dat.res.M0;

Brightness = get(h.slider_Sat,'Value');
h.dat.graph.img1.V =BG+Brightness*h.dat.res.M0.*~h.dat.graph.img0.BackGround...
             .* h.dat.cl.k1;
h.dat.graph.img2.V =BG+Brightness*h.dat.res.M0.*~h.dat.graph.img0.BackGround...
             .*~h.dat.cl.k1;
       
switch Tag
   case {'tg_DisplayMode','tg_EllipseMode','tg_ROIMergeMode','tg_ROISplitMode','tb_LabelMode'}
    
      % do nothing more. 
               
    case 'tb_SVMMode'
        
        for ii=1:h.dat.ops.Nk
            %             h.dat.cellmask(h.dat.res.iclust==h.dat.F.ichosen_append(ii))=1;
            edge_pix = h.dat.stat(ii).ipix_edge;
            type_index = find([h.dat.cl.type_number_table{:,1}]==h.dat.cl.svm_class(ii));
            HSV=rgb2hsv(h.dat.cl.type_number_table{type_index,3});
            
            if type_index ~=1
                  h.dat.graph.img1.V(edge_pix)=HSV(3);
            end
          
%             h.dat.graph.img2.V(edge_pix)=HSV(3);
        end
  
        
    otherwise
        error('Unknown Tag %s',Tag);
end
    
% iselect = h.dat.res.iclust==h.dat.F.ichosen;
% h.dat.graph.img1.V(iselect) = h.dat.cl.k1(iselect);
% h.dat.graph.img2.V(iselect) = ~h.dat.cl.k1(iselect);


for ii=h.dat.F.ichosen_append(:)'
    edge_pix = h.dat.stat(ii).ipix_edge;
   
    h.dat.graph.img1.V(edge_pix)=1;
    h.dat.graph.img2.V(edge_pix)=1;
end



function Smooth = my_conv_local(S1, sig)
NN = size(S1,1);
NT = size(S1,2);
dt = -4*sig:1:4*sig;
gaus = exp( - dt.^2/(2*sig^2));
gaus = gaus'/sum(gaus);
Smooth = filter(gaus, 1, [S1' ones(NT,1); zeros(4*sig, NN+1)]);
Smooth = Smooth(1+4*sig:end, :);
Smooth = Smooth(:,1:NN) ./ (Smooth(:, NN+1) * ones(1,NN));
Smooth = Smooth';


% --- Executes on key press with focus on pb_rois_display and none of its controls.
function pb_rois_display_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to pb_rois_display (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pb_reset_originalselection.
function pb_reset_originalselection_Callback(hObject, eventdata, handles)
% hObject    handle to pb_reset_originalselection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pb_test.
function pb_test_Callback(hObject, eventdata, handles)
% hObject    handle to pb_test (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
show_corr_roi(handles);
% handles=init_movie_timer(hObject,eventdata,handles);
% guidata(hObject,handles);

% --- Executes on button press in pb_movie_display.
function pb_movie_display_Callback(hObject, eventdata, h)
% hObject    handle to pb_movie_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% h.dat.graph.map = 5;
% h=redraw_meanimg(h); 
h=redraw_movie(hObject,eventdata,h);
guidata(hObject,h);




% --- Executes on button press in pb_loadmovie.
function handles=pb_loadmovie_Callback(hObject, eventdata, handles)
% hObject    handle to pb_loadmovie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function handles=local_loadmovie(hObject,eventdata,handles)
[loadpath,loadname,~]=fileparts(handles.dat.filename);

RegFastFilePath=fullfile(loadpath,sprintf('Plane%d_ch%d',handles.dat.ops.PlaneID,handles.dat.ops.ChannelID),...
    sprintf('x%dmovie',handles.dat.ops.RegFile_xtime));
RegFastFile=sprintf('%s_%s_%s_2P_plane%d_ch%d_x%d.tif',...
    handles.dat.ops.date, handles.dat.ops.SubDirs{1},handles.dat.ops.mouse_name,handles.dat.ops.PlaneID,...
    handles.dat.ops.ChannelID,handles.dat.ops.RegFile_xtime);
RegFastFileFindString=sprintf('%s_%s_%s_2P_plane%d_ch%d_x%d*.tif',...
    handles.dat.ops.date, '*',handles.dat.ops.mouse_name,handles.dat.ops.PlaneID,...
    handles.dat.ops.ChannelID,handles.dat.ops.RegFile_xtime);

handles.dat.reg_decimation = handles.dat.ops.RegFile_xtime;

if exist(fullfile(RegFastFilePath,RegFastFile),'file')
    fprintf('reg-fast file found (%s).\n',fullfile(RegFastFilePath,RegFastFile));
else
    f=dir(fullfile(RegFastFilePath,RegFastFileFindString));
    RegFastFile={f.name};
    if isempty(f)
        error('reg-fast file nor similar files are not found: %s',fullfile(RegFastFilePath,RegFastFileFindString));
    end
end

firstIdx = 1;
lastIdx = [];
stride = 1;
ind_i = 1:max(handles.dat.QV.ylim);
ind_j = 1:max(handles.dat.QV.xlim);

handles.dat.reg_fast_filename=fullfile(RegFastFilePath,RegFastFile);

% from RoiGui_007, use mlttiff class to access multiple Tiff files as a single object. 
    handles.dat.reg_data=mlttiff(handles.dat.reg_fast_filename);
% 
% handles.dat.reg_data = loadFramesBuff2(handles.dat.reg_fast_filename, ...
%      firstIdx, lastIdx, stride, ...
%     ind_i,ind_j);

% handles.mij=mij.createImage(handles.dat.reg_fast_filename);
% pathname=strrep(sprintf('path=[%s]',handles.dat.reg_fast_filename),'\','\\')

% memory error in this code... 
% handles.mij=MIJ.run('Open...', pathname);
guidata(hObject,handles);


% --- Executes on button press in tg_DisplayMode.
function handles=tg_DisplayMode_Callback(hObject, eventdata, handles)
% hObject    handle to tg_DisplayMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tg_DisplayMode

% --- Executes on button press in tg_EllipseMode.
function handles=tg_EllipseMode_Callback(hObject, eventdata, handles)
% hObject    handle to tg_EllipseMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tg_EllipseMode


function TimeSelectionButtonDownFcn(hObject,eventdata,h)
h=guidata(hObject);
pos=get(gca,'CurrentPoint');

% eventdata
% z = round(eventdata.Source.CurrentAxes.CurrentPoint(1,:));
frame_index = round(pos(1,1));
y  = round(pos(1,2));

SelectionType =get(gcf,'SelectionType');
% display(SelectionType)

switch SelectionType
    case 'alt' % Right click
    case 'normal' % left click
            fprintf('Frame=%d\n',frame_index);
            goto_frame(hObject,eventdata,h,frame_index)
            
    case 'open' % double click
     
    case 'extend' % shift left
     
end



% --- Executes on mouse press over axes background.
function axes_fluorescence_ButtonDownFcn(hObject, eventdata, h)
% hObject    handle to axes_fluorescence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
TimeSelectionButtonDownFcn(hObject,eventdata,h);



function edit_minSkewF_Callback(hObject, eventdata, h)
% hObject    handle to edit_minSkewF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_minSkewF as text
%        str2double(get(hObject,'String')) returns contents of edit_minSkewF as a double

h.dat.cl.skewF_low = str2double(get(h.edit_minSkewF,'String'));
h = ApplyROIFilter(h);
h = buildLambdaValue(h);
h=redraw_figure(h);
guidata(hObject,h);


% --- Executes during object creation, after setting all properties.
function edit_minSkewF_CreateFcn(hObject, eventdata, h)
% hObject    handle to edit_minSkewF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --------------------------------------------------------------------
function handles=ui_imellipse_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to ui_imellipse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles.dat,'ImEllipseH')
    n = length(handles.dat.ImEllipseH);
else
    n = 0;
end
handles.dat.ImEllipseH(n+1)=imellipse(handles.axes_left);
guidata(hObject,handles);

function h=add_imellipse_Callback(hObject, eventdata, h)
% hObject    handle to ui_imellipse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(h.dat,'ImEllipseH')
    n = length(h.dat.ImEllipseH);
else
    n = 0;
end
h.dat.ImEllipseH(n+1)=imellipse(h.axes_left);
% waitfor(h.dat.ImEllipseH(n+1));

guidata(hObject,h);



% --- Executes on button press in pb_show_EllipticROI.
function handles=pb_show_EllipticROI_Callback(hObject, eventdata, handles)
% hObject    handle to pb_show_EllipticROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

IEh=handles.dat.ImEllipseH;
if isfield(handles.dat,'reg_data')
    % do nothing
else
     handles=local_loadmovie(hObject, eventdata, handles);
     handles.dat.reg_data.get_single_slice(1);
end
  [H,W,D]=size(handles.dat.reg_data.buffer);
  
  handles.dat.ImEllipse=[];
cnt=1;
data=reshape(handles.dat.reg_data.buffer,H*W,[]);
for ii=1:length(IEh)
    if isobject(IEh(ii))
        
%         myfigure('Check mask');
%         imagesc(IEh(ii).createMask);
        [ind_y,ind_x]=find(IEh(ii).createMask);
        ind_y=handles.dat.ops.yrange(ind_y);
        ind_x=handles.dat.ops.xrange(ind_x);
        EllipROI_index=sub2ind([H,W],ind_y,ind_x);
%         index = bsxfun(@plus,EllipROI_index,[0:W*H:W*H*(D-1)]);
        handles.dat.ImEllipse(cnt).EllipROI_index = EllipROI_index;
        handles.dat.ImEllipse(cnt).F_Ellip = mean(data(EllipROI_index,:),1)';
        cnt=cnt+1;
    end
end

myfigure('ImEllipse F');clf;
plot(1:D,[handles.dat.ImEllipse(:).F_Ellip]);

guidata(hObject,handles);

% --- Executes on button press in pb_DeleteEllipticROIs.
function pb_DeleteEllipticROIs_Callback(hObject, eventdata, handles)
% hObject    handle to pb_DeleteEllipticROIs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if isfield(handles.dat,'ImEllipseH')
   delete(handles.dat.ImEllipseH)
   handles.dat=rmfield(handles.dat,'ImEllipseH');
end
guidata(hObject,handles);


% --- Executes on button press in pb_CorrROI_display.
function h=pb_CorrROI_display_Callback(hObject, eventdata, h)
% hObject    handle to pb_CorrROI_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h.dat.graph.map = 1;
h=redraw_corrimg(h);
guidata(hObject,h);


% --- Executes during object creation, after setting all properties.
function pb_CorrROI_display_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pb_CorrROI_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on scroll wheel click while the figure is in focus.
function figure1_WindowScrollWheelFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	VerticalScrollCount: signed integer indicating direction and number of clicks
%	VerticalScrollAmount: number of lines scrolled for each click
% handles    structure with handles and user data (see GUIDATA)

[eventdata.VerticalScrollCount, eventdata.VerticalScrollAmount]


% --------------------------------------------------------------------
function uipanel46_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to uipanel46 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pb_ROIMergeMode.
function pb_ROIMergeMode_Callback(hObject, eventdata, handles)
% hObject    handle to pb_ROIMergeMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pb_ROIMergeMode
if  get(hObject,'Value')
    merge_rois(hObject,eventdata,h)
end


function merge_rois(hObject,eventdata,h)
ROIs_to_merge=h.dat.F.ichosen_append

% --- Executes on button press in tg_ROISplitMode.
function tg_ROISplitMode_Callback(hObject, eventdata, handles)
% hObject    handle to tg_ROISplitMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tg_ROISplitMode


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over tg_DisplayMode.
function tg_DisplayMode_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to tg_DisplayMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when selected object is changed in ModeSelectionButtonGroup.
function handles=ModeSelectionButtonGroup_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in ModeSelectionButtonGroup 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Tag=get(hObject,'Tag');
fprintf('%s\n',Tag);

try delete(h.svm_editnow), catch, end 
 
% set(handles.tb_SVMMode,'CData',);
switch Tag
    case 'tg_DisplayMode'
            
        ButtonDownFcn=@(hObject,eventdata)RoiGui_007('CellSelectionButtonDownFcn',...
            hObject,eventdata,guidata(hObject));
        set(handles.left_imageH,'ButtonDownFcn',ButtonDownFcn);
         set(handles.uipanel_Left,'Title','DisplayMode');
    case 'tg_EllipseMode'
        
        ButtonDownFcn=@(hObject,eventdata)RoiGui_007('EllipseSelectionButtonDownFcn',...
            hObject,eventdata,guidata(hObject));
        set(handles.left_imageH,'ButtonDownFcn',ButtonDownFcn);
        set(handles.uipanel_Left,'Title','EllipseMode');
%     case 'tg_ROIMergeMode' % merge is not considered as a mode. Merge already selected rois. 
       
        
    case 'tg_ROISplitMode'
        
        ButtonDownFcn=@(hObject,eventdata)RoiGui_007('CellSplitButtonDownFcn',...
            hObject,eventdata,guidata(hObject));
        set(handles.left_imageH,'ButtonDownFcn',ButtonDownFcn);
        set(handles.uipanel_Left,'Title','ROISplitMode');
    case 'tb_SVMMode'
        
        ButtonDownFcn=@(hObject,eventdata)RoiGui_007('CellSelectionButtonDownFcn',...
            hObject,eventdata,guidata(hObject));
%         ButtonDownFcn=@(hObject,eventdata)RoiGui_007('SVMSelectionButtonDownFcn',...
%             hObject,eventdata,guidata(hObject));
        set(handles.left_imageH,'ButtonDownFcn',ButtonDownFcn);
         set(handles.uipanel_Left,'Title','SVMMode');
    case 'tb_LabelMode'
         ButtonDownFcn=@(hObject,eventdata)RoiGui_007('CellSelectionButtonDownFcn',...
            hObject,eventdata,guidata(hObject));
        set(handles.left_imageH,'ButtonDownFcn',ButtonDownFcn);
         set(handles.uipanel_Left,'Title','DisplayMode');
    otherwise
        error('Unknown Tag %s',Tag);
end
update_figure(handles);

guidata(hObject,handles);

function h=renew_cluster(h)
h.dat.cl.IDX=recalc_IDX(h);
h=update_iclust1(h);

function  IDX=recalc_IDX(h)
epsilon = str2double(  get(h.edit_CorrROI_epsilon ,'String'));
MinPts = str2double(  get(h.edit_CorrROI_MinPts,'String'));
IDX=dbscan_D(1-h.dat.cl.C_of_zF,epsilon,MinPts);

unique_cluster=unique(IDX);
unique_cluster(unique_cluster==0)=[];

IDXtmp = zeros(size(IDX));

for ii=1:length(unique_cluster)
    ind = find(IDX==unique_cluster(ii));
    npix=h.dat.cl.statTBL.npix(ind);
    if iscell(npix)
        npix=[npix{:}];
    end
    [npix_max,npix_max_ind]=  max(npix);
    IDXtmp(ind)= ind(npix_max_ind);
end
IDX = IDXtmp;

function h=update_iclust1(h)

% IDX=h.dat.cl.IDX;
% h.dat.res.iclust1 = h.dat.res.iclust;
% 
% % first, assign index to each cluster
% h.dat.cl.automerge =IDX;
% ind = find(h.dat.cl.automerge==0);
% h.dat.cl.automerge(ind)=ind; 
% 
% % overwrite the cluster
% ind=find(h.dat.cl.manualmerge);
% h.dat.cl.roi_cluster = h.dat.cl.automerge;
% h.dat.cl.roi_cluster(ind)=h.dat.cl.manualmerge(ind);
% 
% unique_cluster=unique(h.dat.cl.roi_cluster);
% unique_cluster(unique_cluster==0)=[];
% 
% for ii=1:length(unique_cluster)
%     ind = find(IDX==unique_cluster(ii));
%     [npix_max,npix_max_ind]=  max(h.dat.cl.statTBL.npix(ind));
%     for jj=1:length(ind)
%     h.dat.res.iclust1(h.dat.res.iclust1==ind(jj))=ind(npix_max_ind); % set the ID to the ROI with biggest npix
%     end
% end

% --- Executes on selection change in listbox_ROIs.
function h=listbox_ROIs_Callback(hObject, eventdata, h)
% hObject    handle to listbox_ROIs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_ROIs contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_ROIs
contents = cellstr(get(hObject,'String'));
contents{get(hObject,'Value')};
ROIindex = get(hObject,'UserData');
selected_ROI_index = ROIindex(get(hObject,'Value'));

Tag=h.ModeSelectionButtonGroup.SelectedObject.Tag;


switch Tag
    case 'tg_DisplayMode'
     % do nothing
        
    case 'tg_EllipseMode'
    % do nothing
       
    case 'tg_ROIMergeMode'
     % convert cluster to roi elements index  
     all_candidates = find(h.dat.cl.automerge==selected_ROI_index);
     included_candidates =find(h.dat.cl.roi_cluster ==selected_ROI_index);
     
     excluded_candidate = setdiff(all_candidates,included_candidates);
    
     selected_ROI_index= included_candidates(1);
     val=max(1,get(h.listbox_automerge,'Value'));
     set(h.listbox_automerge,'String',num2str(included_candidates),'Max',length(included_candidates),...
         'Value',min(val,length(included_candidates))); 
 
     val=max(1,get(h.listbox_exclude_merge,'Value'));
     set(h.listbox_exclude_merge,'String',num2str(excluded_candidate),'Max',length(excluded_candidate),...
         'Value',min(val,length(excluded_candidate))); 
    otherwise
        error('Unknown Tag %s',Tag);
end

h.dat.F.ichosen = selected_ROI_index;
 h.dat.F.ichosen_append =h.dat.F.ichosen; % reset
 
h=update_figure(h);
guidata(hObject,h);

% --- Executes during object creation, after setting all properties.
function listbox_ROIs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_ROIs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in pb_DetectCorrROIs.
function pb_DetectCorrROIs_Callback(hObject, eventdata, h)
% hObject    handle to pb_DetectCorrROIs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Run DBSCAN Clustering Algorithm

epsilon = str2double(get(h.edit_CorrROI_epsilon,'String'));
MinPts  = str2double(get(h.edit_CorrROI_MinPts, 'String'));

% epsilon=0.5;
% MinPts=2;
D = 1-h.dat.cl.C_of_zF;
% compute cluster index IDX
h.dat.cl.IDX=recalc_IDX(h);
 
unique_cluster = unique(h.dat.cl.IDX);
unique_cluster(unique_cluster==0)=[];
unique_cluster= unique_cluster(:)';
%% Plot Results
% 
% rands_tmp = h.dat.cl.rands;
% 
% for ii=1:length(unique_cluster)
%     index=find(IDX==unique_cluster(ii));
%   for jj=1:length(index)
%       rands_tmp(index)= rands_tmp(index(1));
%   end 
% end
% 
%     
% Smap = ones(size(h.dat.graph.img1.Sat)); % saturation, 1 is more pure color, 0 is white.
% % Hmap = zeros(1,length(h.dat.res.iclust)); % hue color map, make the same cell same color.
% % Hmap(roi_samecell)=h.dat.cl.rands(h.dat.F.ichosen);
% Hmap       = h.dat.graph.img2.H; reshape(Hmap(h.dat.res.iclust), h.dat.res.Ly, h.dat.res.Lx);
% 
% Vmap = h.dat.graph.img2.V; % value, 1 is brighter
% 
% 
% I = hsv2rgb(cat(3, Hmap, Smap, Vmap));
% I = min(I, 1);
% image(I);
Vmap = h.dat.graph.img1.V;
figh=myfigure('CorrROI');clf;% set(figh,'Position',);
for ii=unique_cluster
    
    subplot(2,1,1);
    roi_samecell = find(h.dat.cl.IDX==ii);
    opt.hue_range= [0.3 0.31];
    opt.highlight_labels = roi_samecell; % a vector of indices. ROI belongs to this indices are colored with opt.highlight_color. 
    opt.highlight_color = ones(1,length(roi_samecell));
    
    I=ROI_gem_img(h.dat.res.iclust1,...
                        h.dat.res.probabilities,...
                        h.dat.res.M0,...
                        opt);
    imagesc(I);
    title_txt = sprintf('Cluster=%d ',ii);
    title(title_txt)
    subplot(2,1,2);
    
    y_trace=zscore(my_conv_local(double(h.dat.F.trace(roi_samecell,:)), 3),0,2);
    y_trace = bsxfun(@plus, y_trace,2*[0:length(roi_samecell)-1]');
  
    plot(1:size(y_trace,2),y_trace);
    if ishandle(figh)
        pause
    else
        break;
    end
end

function edit_CorrROI_epsilon_Callback(hObject, eventdata, h)
% hObject    handle to edit_CorrROI_epsilon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_CorrROI_epsilon as text
%        str2double(get(hObject,'String')) returns contents of edit_CorrROI_epsilon as a double

h.dat.cl.IDX=recalc_IDX(h);
guidata(hObject,h);

% --- Executes during object creation, after setting all properties.
function edit_CorrROI_epsilon_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_CorrROI_epsilon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_CorrROI_MinPts_Callback(hObject, eventdata, h)
% hObject    handle to edit_CorrROI_MinPts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_CorrROI_MinPts as text
%        str2double(get(hObject,'String')) returns contents of edit_CorrROI_MinPts as a double
h.dat.cl.IDX=recalc_IDX(h);
guidata(hObject,h);

% --- Executes during object creation, after setting all properties.
function edit_CorrROI_MinPts_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_CorrROI_MinPts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_automerge.
function listbox_automerge_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_automerge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_automerge contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_automerge
contents = cellstr(get(hObject,'String'));
contents{get(hObject,'Value')};

% --- Executes during object creation, after setting all properties.
function listbox_automerge_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_automerge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_exclude_merge.
function listbox_exclude_merge_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_exclude_merge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_exclude_merge contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_exclude_merge


% --- Executes during object creation, after setting all properties.
function listbox_exclude_merge_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_exclude_merge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_ROIMergeSelection.
function h=popup_ROIMergeSelection_Callback(hObject, eventdata, h)
% hObject    handle to popup_ROIMergeSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_ROIMergeSelection contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_ROIMergeSelection
contents = cellstr(get(hObject,'String'));
Mode=contents{get(hObject,'Value')};

switch Mode
    case 'All elements'
        unique_ROI_elements = unique(h.dat.res.iclust1);     
        s=cell(length(unique_ROI_elements),1);
        for ii=1:length(unique_ROI_elements)
            if h.dat.cl.selected(ii)
                s{ii}=sprintf('<HTML><BODY text=%s>%d','black',unique_ROI_elements(ii));
            else
                s{ii}=sprintf('<HTML><BODY text=%s>%d','aaaaaa',unique_ROI_elements(ii));
            end
        end
        set(h.listbox_ROIs,'String',s);
        set(h.listbox_ROIs,'UserData',unique_ROI_elements);
%         get(h.listbox_ROIs)
    case 'merge_candidate'
        unique_cluster=unique(h.dat.cl.IDX);
        unique_cluster(unique_cluster==0)=[];
        s=cell(length(unique_cluster),1);
        for ii=1:length(unique_cluster)
            n=nnz(h.dat.cl.IDX==unique_cluster(ii));
            s{ii}=sprintf('<HTML><BODY text=%s>%d (%d)','black',unique_cluster(ii),n);
        end
        set(h.listbox_ROIs,'String',s,'Value',1 );
        set(h.listbox_ROIs,'UserData',unique_cluster);
        
    case 'merged'
        
end

% --- Executes during object creation, after setting all properties.
function popup_ROIMergeSelection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_ROIMergeSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_exclude_mergedroi.
function pb_exclude_mergedroi_Callback(hObject, eventdata, h)
% hObject    handle to pb_exclude_mergedroi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Str=get(h.listbox_automerge,'String');
ExcludedROIs = str2num(Str(get(h.listbox_automerge,'Value'),:));
h.dat.cl.manualmerge(ExcludedROIs)=ExcludedROIs;
h=update_iclust1(h);

h=listbox_ROIs_Callback(h.listbox_ROIs, eventdata, h); % guidata updated in there

% guidata(hObject,h);

% --- Executes on button press in pb_include_mergedroi.
function pb_include_mergedroi_Callback(hObject, eventdata, h)
% hObject    handle to pb_include_mergedroi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Str=get(h.listbox_exclude_merge,'String');
ToBeIncludedROIs = str2num(Str(get(h.listbox_exclude_merge,'Value'),:));
h.dat.cl.manualmerge(ToBeIncludedROIs)=0;
h=update_iclust1(h);

h=listbox_ROIs_Callback(h.listbox_ROIs, eventdata, h); % guidata updated in there


% --- Executes on button press in rb_iscell_filter.
function rb_iscell_filter_Callback(hObject, eventdata, handles)
% hObject    handle to rb_iscell_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rb_iscell_filter


% --- Executes on slider movement.
function slider_BGLevel_Callback(hObject, eventdata, handles)
% hObject    handle to slider_BGLevel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles = buildLambdaValue(handles); %
handles=redraw_figure(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function slider_BGLevel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_BGLevel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_Sat_Callback(hObject, eventdata, handles)
% hObject    handle to slider_Sat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles = buildLambdaValue(handles); %
handles=redraw_figure(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function slider_Sat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_Sat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end




function edit_EccentricityMax_Callback(hObject, eventdata, h)
% hObject    handle to edit_EccentricityMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_EccentricityMax as text
%        str2double(get(hObject,'String')) returns contents of edit_EccentricityMax as a double

h.dat.cl.EccentricityMax = str2double(get(h.edit_EccentricityMax,'String'));
h = ApplyROIFilter(h);
h = buildLambdaValue(h);
h=redraw_figure(h);
guidata(hObject,h);

% --- Executes during object creation, after setting all properties.
function edit_EccentricityMax_CreateFcn(hObject, eventdata, h)
% hObject    handle to edit_EccentricityMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_SolidityMin_Callback(hObject, eventdata, h)
% hObject    handle to edit_SolidityMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_SolidityMin as text
%        str2double(get(hObject,'String')) returns contents of edit_SolidityMin as a double


h.dat.cl.statTBL.SolidityMin = str2double(get(h.edit_SolidityMin,'String'));
h = ApplyROIFilter(h);
h = buildLambdaValue(h);
h=redraw_figure(h);
guidata(hObject,h);


% --- Executes during object creation, after setting all properties.
function edit_SolidityMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_SolidityMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function h=update_svm_class(hObject,eventdata,h,varargin)

if nargin>=4
    roi_type = varargin{1};
else
    error();
end
% roi_type
% h.dat.F.ichosen_append

h.dat.cl.svm_class(h.dat.F.ichosen_append)=roi_type;
guidata(hObject,h);

% --- Executes on button press in pb_SVM_1_CellBody.
function pb_SVM_1_CellBody_Callback(hObject, eventdata, handles)
% hObject    handle to pb_SVM_1_CellBody (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

update_svm_class(hObject,eventdata,handles,1);

% --- Executes on button press in pb_2_Dendrite.
function pb_2_Dendrite_Callback(hObject, eventdata, handles)
% hObject    handle to pb_2_Dendrite (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
update_svm_class(hObject,eventdata,handles,2);

% --- Executes on button press in pb_3_Axon.
function pb_3_Axon_Callback(hObject, eventdata, handles)
% hObject    handle to pb_3_Axon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
update_svm_class(hObject,eventdata,handles,3);

% --- Executes on button press in pb_SVM_4_Spine.
function pb_SVM_4_Spine_Callback(hObject, eventdata, handles)
% hObject    handle to pb_SVM_4_Spine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
update_svm_class(hObject,eventdata,handles,4);

% --- Executes on button press in pb_SVM_0_Unknown.
function pb_SVM_0_Unknown_Callback(hObject, eventdata, handles)
% hObject    handle to pb_SVM_0_Unknown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
update_svm_class(hObject,eventdata,handles,0);


% --- Executes on button press in pb_SVM_9_Noise.
function pb_SVM_9_Noise_Callback(hObject, eventdata, handles)
% hObject    handle to pb_SVM_9_Noise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

update_svm_class(hObject,eventdata,handles,9);

% --- Executes on button press in pb_ShowSVMTeacher.
function pb_ShowSVMTeacher_Callback(hObject, eventdata, handles)
% hObject    handle to pb_ShowSVMTeacher (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pb_UpdateSVMPrediction.
function pb_UpdateSVMPrediction_Callback(hObject, eventdata, h)
% hObject    handle to pb_UpdateSVMPrediction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

labels = h.dat.cl.svm_class;

% labels = labels(teacher_ind);

svm_axis = {'Compactness','npix','Eccentricity','Solidity','V','skewF','SNratio'};
data = h.dat.cl.statTBL(:,svm_axis);

if isfield(h.dat.cl,'svm_teacher_label')
    loaded_labels = h.dat.cl.svm_teacher_label;
    labels = cat(1,labels,loaded_labels);
else
    loaded_labels= [];
end

unique_classes = unique(labels); unique_classes(unique_classes==0)=[];
teacher_ind = (labels~=0);
if isfield(h.dat.cl,'svm_teacher_data')
    loaded_teacher_data = h.dat.cl.svm_teacher_data;
    data = cat(1,data,loaded_teacher_data);
else
    loaded_teacher_data= [];
end
N_imported_data = length(loaded_labels);
data = nanzscore(double(table2array(data)),0,1);

% data= double(data(:,teacher_ind));

fprintf('training...\n')

for cc=unique_classes(:)'
    train_label= double(labels==cc);
    
    fprintf('class=%d:(N=%d)\n',cc,nnz(train_label));
    h.dat.cl.svmmodel(cc).param = svmtrain(2*train_label(teacher_ind)-1,...
        data(teacher_ind,:), '-c 400 -g 0.01');
    
    [predict_label, accuracy, dec_values] = svmpredict(randn(size(data,1),1), ...
        data,  h.dat.cl.svmmodel(cc).param); % test the training data
    
    h.dat.cl.svm_predicted_type(predict_label>0)=cc;
end
% overwrite with the teacher data
local_teacher_ind = teacher_ind(1:end-length(loaded_labels));
% h.dat.cl.svm_predicted_type(local_teacher_ind)=h.dat.cl.svm_class(local_teacher_ind);
 h.dat.cl.svm_predicted_type(local_teacher_ind)=h.dat.cl.svm_class(local_teacher_ind);
 h.dat.cl.svm_predicted_type = h.dat.cl.svm_predicted_type(1:end-length(loaded_labels));
 h = ApplyROIFilter(h);
h = buildHue(h);
h = buildSat(h);
h = buildLambdaValue(h);

h=update_figure(h);            
guidata(hObject,h);
        

% --- Executes on button press in rb_UseSVMPredictionAsFilter.
function rb_UseSVMPredictionAsFilter_Callback(hObject, eventdata, handles)
% hObject    handle to rb_UseSVMPredictionAsFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rb_UseSVMPredictionAsFilter


% --- Executes on button press in cb_SVM_CellBody.
function cb_SVM_CellBody_Callback(hObject, eventdata, handles)
% hObject    handle to cb_SVM_CellBody (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_SVM_CellBody


% --- Executes on button press in cb_SVM_Dendrite.
function cb_SVM_Dendrite_Callback(hObject, eventdata, handles)
% hObject    handle to cb_SVM_Dendrite (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_SVM_Dendrite


% --- Executes on button press in cb_SVM_Axon.
function cb_SVM_Axon_Callback(hObject, eventdata, handles)
% hObject    handle to cb_SVM_Axon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_SVM_Axon


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4


% --- Executes on selection change in popup_FilterSelection.
function popup_FilterSelection_Callback(hObject, eventdata, handles)
% hObject    handle to popup_FilterSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_FilterSelection contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_FilterSelection


% --- Executes during object creation, after setting all properties.
function popup_FilterSelection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_FilterSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in rb_BorderTeacherROI.
function rb_BorderTeacherROI_Callback(hObject, eventdata, handles)
% hObject    handle to rb_BorderTeacherROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rb_BorderTeacherROI


% --- Executes on button press in pb_loadSVM.
function pb_loadSVM_Callback(hObject, eventdata, h)
% hObject    handle to pb_loadSVM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(h, 'dat') && isfield(h.dat, 'filename')
    root = fileparts(h.dat.filename);
else
    root = 'G:Kosuk\DATA\F\';
end
[filename1,filepath1]=uigetfile(root, 'Select _PROC.mat File to load SVM parameters');
LOADNAME = fullfile(filepath1, filename1);
fprintf('Loading %s...\n',LOADNAME);
load(LOADNAME,'cl');

labels = cl.svm_class;
teacher_ind = (labels~=0);
% labels = labels(teacher_ind);
unique_classes = unique(labels); unique_classes(unique_classes==0)=[];

svm_axis = {'Compactness','npix','Eccentricity','Solidity','V','skewF','SNratio'};
teacher_data = cl.statTBL(teacher_ind,svm_axis);
h.dat.cl.svm_teacher_label = labels(teacher_ind);
h.dat.cl.svm_teacher_data = teacher_data;
% teacher_data = 
% h.dat.cl.svm_class = cl.svm_class; 
guidata(hObject,h);


% --- Executes on button press in pb_Print.
function pb_Print_Callback(hObject, eventdata, h)
% hObject    handle to pb_Print (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

printname = h.dat.filename(1:end-4);
fprintf('Printing %s ...\n',printname);
set(gcf,'PaperPositionMode','auto');
print(gcf,'-dpng',[printname,'.png'])
print(gcf,'-depsc2',[printname,'.eps']);


% --- Executes on button press in pb_ResetSVM.
function pb_ResetSVM_Callback(hObject, eventdata, h)
% hObject    handle to pb_ResetSVM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


h.dat.cl.svm_class = zeros(size(h.dat.cl.svm_class));

labels = h.dat.cl.svm_class;

teacher_ind = [];
% labels = labels(teacher_ind);
unique_classes = [];

svm_axis = {'Compactness','npix','Eccentricity','Solidity','V','skewF','SNratio'};
data = h.dat.cl.statTBL(:,svm_axis);

if isfield(h.dat.cl,'svm_teacher_label')
labels = cat(1,labels,h.dat.cl.svm_teacher_label);
end

if isfield(h.dat.cl,'svm_teacher_data')
data = cat(1,data,h.dat.cl.svm_teacher_data);
end
data = zscore(double(table2array(data)),0,1);


% data= double(data(:,teacher_ind));

fprintf('training...\n')

h.dat.cl.svm_predicted_type= zeros(size(labels));
h.dat.cl.manual = zeros(size(h.dat.cl.manual));

for cc=unique_classes(:)'
    train_label= double(labels==cc);
    
    fprintf('class=%d:(N=%d)\n',cc,nnz(train_label));
    h.dat.cl.svmmodel(cc).param = svmtrain(2*train_label(teacher_ind)-1,...
        data(teacher_ind,:), '-c 400 -g 0.01');
    
    [predict_label, accuracy, dec_values] = svmpredict(randn(h.dat.ops.Nk,1), ...
        data,  h.dat.cl.svmmodel(cc).param); % test the training data
    
    h.dat.cl.svm_predicted_type(predict_label>0)=cc;
end
% overwrite with the teacher data
h.dat.cl.svm_predicted_type(teacher_ind)=h.dat.cl.svm_class(teacher_ind);
h = ApplyROIFilter(h);
h = buildHue(h);
h = buildSat(h);
h = buildLambdaValue(h);

h=update_figure(h);            
guidata(hObject,h);



% --- Executes on selection change in popup_ROI_ColorSelection.
function popup_ROI_ColorSelection_Callback(hObject, eventdata, handles)
% hObject    handle to popup_ROI_ColorSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_ROI_ColorSelection contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_ROI_ColorSelection


% --- Executes during object creation, after setting all properties.
function popup_ROI_ColorSelection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_ROI_ColorSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cb_ROI_show.
function cb_ROI_show_Callback(hObject, eventdata, handles)
% hObject    handle to cb_ROI_show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_ROI_show


% --- Executes on selection change in popup_Label1_ColorSelection.
function popup_Label1_ColorSelection_Callback(hObject, eventdata, handles)
% hObject    handle to popup_Label1_ColorSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_Label1_ColorSelection contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_Label1_ColorSelection


% --- Executes during object creation, after setting all properties.
function popup_Label1_ColorSelection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_Label1_ColorSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cb_show_Label1.
function cb_show_Label1_Callback(hObject, eventdata, handles)
% hObject    handle to cb_show_Label1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_show_Label1


% --- Executes on selection change in popup_Label2_ColorSelection.
function popup_Label2_ColorSelection_Callback(hObject, eventdata, handles)
% hObject    handle to popup_Label2_ColorSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_Label2_ColorSelection contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_Label2_ColorSelection


% --- Executes during object creation, after setting all properties.
function popup_Label2_ColorSelection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_Label2_ColorSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cb_show_Label2.
function cb_show_Label2_Callback(hObject, eventdata, handles)
% hObject    handle to cb_show_Label2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_show_Label2


% --- Executes on button press in checkbox9.
function checkbox9_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox9


% --- Executes on button press in checkbox10.
function checkbox10_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox10


% --- Executes on button press in pb_LoadLabelImage.
function pb_LoadLabelImage_Callback(hObject, eventdata, h)
% hObject    handle to pb_LoadLabelImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

RegImg= Load_register_OtherChannels_20181001(h);

% --- Executes on button press in tb_LabelMode.
function tb_LabelMode_Callback(hObject, eventdata, handles)
% hObject    handle to tb_LabelMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tb_LabelMode


% --- Executes during object creation, after setting all properties.
function rb_showneuropil_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rb_showneuropil (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
