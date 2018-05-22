function varargout = RoiGui_004(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RoiGui_004_OpeningFcn, ...
                   'gui_OutputFcn',  @RoiGui_004_OutputFcn, ...
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

% --- Executes just before RoiGui_004 is made visible.
function RoiGui_004_OpeningFcn(hObject, eventdata, h, varargin)
% This function has no output args, see OutputFcn.
% varargin   command line arguments to RoiGui_004 (see VARARGIN)
 
h=init_icons(hObject,eventdata,h);
h=init_movie_timer(hObject,eventdata,h);
h.control_on = 0;
h.shift_on = 0;
 
h.output = hObject;
guidata(hObject, h);
% UIWAIT makes RoiGui_004 wait for user response (see UIRESUME)
% uiwait(h.figure1);

function h=init_icons(hObject,eventdata,h)
h.icons.SVM_BW  =imread('C:\home\matlab_svn\Suite2P\KHcodes\icons\SVM_Black_26x26.png');
h.icons.SVM_Color =imread('C:\home\matlab_svn\Suite2P\KHcodes\icons\SVM_26x26.png');
h.icons.ROICircle_BW=imread('C:\home\matlab_svn\Suite2P\KHcodes\icons\ROICircle_25x26_Gray.png');
h.icons.ROICircle_Color=imread('C:\home\matlab_svn\Suite2P\KHcodes\icons\ROICircle_25x26_Color.png');
h.icons.CellDisplay_BW=imread('C:\home\matlab_svn\Suite2P\KHcodes\icons\CellDisplay_25x26_Gray.png');
h.icons.CellDisplay_Color=imread('C:\home\matlab_svn\Suite2P\KHcodes\icons\CellDisplay_25x26_Color.png');


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
function varargout = RoiGui_004_OutputFcn(hObject, eventdata, h) 
% varargout  cell array for returning output args (see VARARGOUT);
varargout{1} = h.output;

function pb_LoadPlane_Callback(hObject, eventdata, h)
% [filename1,filepath1]=uigetfile('\\zserver\Lab\Share\Marius\', 'Select Data File');

flag = 0;
try
    if isfield(h, 'dat') && isfield(h.dat, 'filename')
        root = fileparts(h.dat.filename);
    else
        root = 'G:Kosuk\DATA\F\';
    end
    [filename1,filepath1]=uigetfile(root, 'Select Data File');
    h.dat = load(fullfile(filepath1, filename1));
    set(h.figure1, 'Name', filename1);

    flag = 1;
catch
end

if flag
    % if the user selected a file, do all the initializations
    rng('default')


% keyboard;
if isfield(h.dat, 'dat')
    h.dat = h.dat.dat;
else
    h.dat.filename = fullfile(filepath1, filename1);
    
    h.dat.cl.Mrs      = [h.dat.stat.mrs]./[h.dat.stat.mrs0]; % stat0 (stat0.mrs = inf) 
    h.dat.cl.npix     = [h.dat.stat.npix];
    h.dat.cl.Ly       = numel(h.dat.ops.yrange);
    h.dat.cl.Lx       = numel(h.dat.ops.xrange);
    h.dat.cl.MeanM    = 2*mean(h.dat.res.M);
    h.dat.cl.excluded_pixels  = zeros(h.dat.cl.Ly, h.dat.cl.Lx);
    h.dat.cl.excluded_regions = zeros(h.dat.cl.Ly, h.dat.cl.Lx);
    h.dat.cl.excl_pix_perc    = zeros(h.dat.cl.Ly, h.dat.cl.Lx);
    h.dat.cl.topregion        = ones(h.dat.cl.Ly, h.dat.cl.Lx);
    h.dat.cl.Eccentricity = [h.dat.stat.Eccentricity];
    h.dat.cl.Solidity = [h.dat.stat.Solidity];
    h.dat.cl.V        = [h.dat.stat.V]; 
    
    if isfield(h.dat.stat, 'parent')
        h = get_parent_stats(h);
    end
    
    h.dat.res.iclust = reshape(h.dat.res.iclust, h.dat.cl.Ly, h.dat.cl.Lx);
%     h.dat.res.BackGround = h.dat.res.iclust==max(h.dat.res.iclust(:));
    h.dat.res.iclust1 = h.dat.res.iclust; % update iclust to merge. 
    
    
    Nk = h.dat.ops.Nk;
    h.dat.ops.Nk = numel(h.dat.stat);
    rng('default');
    h.dat.cl.rands_orig   = .1 + .8 * rand(1, h.dat.ops.Nk);
    h.dat.cl.rands        = h.dat.cl.rands_orig;
    h.dat.cl.manualmerge   = zeros(h.dat.ops.Nk,1);
    h.dat.cl.manual_type   = zeros(h.dat.ops.Nk,1);
    
    % roi_type, 'type', RGB 
    h.dat.cl.type_number_table = ...
        {0,'Unknown',   [1 1 1];...
        1, 'Cell'       [1 0 0]; ...
        2,'Dendrite',   [0 1 0];...
        3, 'Axon',      [0 0 1];...
        4, 'Spine',     [1 1 0];...
        9, 'Noise',     [1 1 1]}; 
     h.dat.cl.predicted_type   = zeros(h.dat.ops.Nk,1);
     
    if isfield(h.dat, 'clustrules')
         % ROI rules
         h.dat.res.Mrs_thresh_orig = h.dat.clustrules.Compact;
         h.dat.cl.npix_low_orig    = h.dat.clustrules.MinNpix;
         h.dat.cl.npix_high_orig   = h.dat.clustrules.MaxNpix;
    else
        % ROI rules
        h.dat.res.Mrs_thresh_orig   = 3;
        h.dat.cl.npix_low_orig      = 20;
        h.dat.cl.npix_high_orig     = 500;
    end

    % parent rules
%     h.dat.cl.mrs_parent_max = Inf;
%     h.dat.cl.npix_res_max   = Inf;
%     h.dat.cl.npix_par_max   = Inf;
%     h.dat.cl.nreg_max       = Inf;
%     h.dat.cl.VperPix_min    = 0;

    h.dat.cl.skewF          =2*ones(1,length(h.dat.cl.npix)); % init skewness 
       
    h = setOriginalThresh(h);
    
    set(h.edit_Compactness,'String', num2str(h.dat.res.Mrs_thresh));
    set(h.edit_PixelCountHigh,'String', num2str(h.dat.cl.npix_high));
    set(h.edit_PixelCountLow,'String', num2str(h.dat.cl.npix_low));
    
    % set all quadrants as not visited
    h.quadvalue = zeros(3);
    for j = 1:3
        for i = 1:3
            set(h.(sprintf('Q%d%d', j,i)), 'BackgroundColor',[.92 .92 .92]);
        end
    end
    % start with unit vector map
    h.dat.display_select = true; % start with selection view
    lam = h.dat.res.lambda;
    h.dat.img0.V = max(0, min(1, .5 * reshape(lam, h.dat.cl.Ly, h.dat.cl.Lx)/mean(lam(:))));
    
    h.dat.ylim = [0 h.dat.cl.Ly];
    h.dat.xlim = [0 h.dat.cl.Lx];    
    
    h.dat.cl.manual  = zeros(h.dat.ops.Nk, 1);
    if ~isfield(h.dat.cl,'redcell')
        h.dat.cl.redcell = zeros(h.dat.ops.Nk, 1);
    end
    
    h                = splitROIleftright(h);
    
    icell = find(h.dat.cl.iscell);
    if ~isempty(icell)
        h.dat.F.ichosen = icell(1); %ceil(rand * numel(icell))
    else
        h.dat.F.ichosen = 1; %ceil(rand * numel(icell))
    end
    h.dat.F.ichosen_append = [];
    
    Sat = ones(h.dat.cl.Ly, h.dat.cl.Lx);
%     Sat(h.dat.res.iclust==h.dat.F.ichosen) = 0;
    h.dat.img0.BackGround = h.dat.res.probabilities;
    Sat(h.dat.img0.BackGround)=0; % make the background gray color.
    h.dat.img1.Sat     = Sat;
    h.dat.img2.Sat     = Sat;
    
    Mmin = min(h.dat.res.M(:));
    Mmax = max(h.dat.res.M(:));
    h.dat.res.M0 = (h.dat.res.M-Mmin+0.05)/(Mmax-Mmin);
        
    % loop through redcells and set h.dat.cl.rands(h.dat.F.ichosen) = 0
    for j = find(h.dat.cl.redcell)
        h.dat.F.ichosen = j;
        h.dat.cl.rands(h.dat.F.ichosen) = 0;
    end
    if ~isempty(icell)
        h.dat.F.ichosen = icell(1); %ceil(rand * numel(icell))
    else
        h.dat.F.ichosen = 1; %ceil(rand * numel(icell))
    end
    h = buildHue(h);
    h = buildLambdaValue(h);
    
    % x and y limits on subquadrants
    h.dat.figure.x0all = round(linspace(0, 19/20*h.dat.cl.Lx, 4));
    h.dat.figure.y0all = round(linspace(0, 19/20*h.dat.cl.Ly, 4));
    h.dat.figure.x1all = round(linspace(1/20 * h.dat.cl.Lx, h.dat.cl.Lx, 4));
    h.dat.figure.y1all = round(linspace(1/20 * h.dat.cl.Ly, h.dat.cl.Ly, 4));
    
    h.dat.F.Fcell = h.dat.Fcell; 
    h.dat.Fcell = [];    
    
    if isfield(h.dat, 'FcellNeu')
        h.dat.F.FcellNeu = h.dat.FcellNeu; h.dat.FcellNeu = [];
        if mean(sign(h.dat.F.FcellNeu{1}(:)))<0
            for j = 1:length(h.dat.F.FcellNeu)
                h.dat.F.FcellNeu{j} = - h.dat.F.FcellNeu{j};
                h.dat.F.Fcell{j} = h.dat.F.Fcell{j} + h.dat.F.FcellNeu{j};
            end
        end    
        
%         if isfield(h.dat.cl, 'dcell')
%             for k = 1:length(h.dat.cl.dcell)
%                 for j = 1:length(h.dat.F.FcellNeu)
%                     if isfield(h.dat.cl.dcell{k}, 'B')
%                         c2 = h.dat.cl.dcell{k}.B(3);
%                         c1 = h.dat.cl.dcell{k}.B(2);
%                         h.dat.F.FcellNeu{j}(k+Nk, :) = c1 + c2 * h.dat.F.FcellNeu{j}(k+Nk, :);
%                     end
%                 end
%             end
%         end
    end
end

h.dat.maxmap = 1;
ops = h.dat.ops;

% h.dat.mimg(:,:,1)  is mean image
if isfield(ops, 'mimg1') && ~isempty(ops.mimg1)
    h.dat.maxmap = h.dat.maxmap + 1;
    h.dat.mimg(:,:,h.dat.maxmap) = ops.mimg1(ops.yrange, ops.xrange);
    h.dat.mimg_proc(:,:,h.dat.maxmap) = normalize_image(h.dat.mimg(:,:,h.dat.maxmap));
end

if isfield(ops, 'mimgRED') && ~isempty(ops.mimgRED)
    h.dat.maxmap = h.dat.maxmap + 1;
    h.dat.mimg(:,:,h.dat.maxmap) = ops.mimgRED(ops.yrange, ops.xrange);
    h.dat.mimg_proc(:,:,h.dat.maxmap) = normalize_image(h.dat.mimg(:,:,h.dat.maxmap));
end

if isfield(ops, 'mimgREDcorrected') && ~isempty(ops.mimgREDcorrected)
    h.dat.maxmap = h.dat.maxmap + 1;
    h.dat.mimg(:,:,h.dat.maxmap) = ops.mimgREDcorrected;
    h.dat.mimg_proc(:,:,h.dat.maxmap) = normalize_image(h.dat.mimg(:,:,h.dat.maxmap));
end

h.dat.procmap = 0;

h=update_Ftrace(h);


% to init image in left/right axes
axes(h.axes_left); 
h.left_imageH=imagesc(h.dat.mimg(:,:,h.dat.map)); 
colormap('gray'); axis off
% set(h.left_imageH,'CDataMapping','direct');

% axes(h.axes_right);  
% h.left_imageH=imagesc(h.dat.mimg(:,:,h.dat.map)); 
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
 
ButtonDownFcn=@(hObject,eventdata)RoiGui_004('TimeSelectionButtonDownFcn',h.axes_fluorescence,eventdata,h);
set(h.axes_fluorescence,'ButtonDownFcn',ButtonDownFcn);

h=renew_cluster(h);
h=redraw_fluorescence(h);
h=redraw_figure(h);

guidata(hObject,h)

end


function obj=init_fluorescence_plot_handles(h, n)
obj = zeros(n,1);

[NN NT] = size(h.dat.F.trace);
for ii=1:n
    obj(ii) = plot(1:NT,zeros(1,NT));
    hold on;
end


function h=update_Ftrace(h)
%% construct F.trace
h.dat.map = 1;
h.dat.F.trace = [];
for i = 1:length(h.dat.F.Fcell)
    h.dat.F.trace = cat(2, h.dat.F.trace, h.dat.F.Fcell{i});
end
if isfield(h.dat.F, 'FcellNeu')
    h.dat.F.neurop = [];
    for i = 1:length(h.dat.F.FcellNeu)
        h.dat.F.neurop = cat(2, h.dat.F.neurop, h.dat.F.FcellNeu{i});
    end    
    
else
   h.dat.F.neurop = zeros(size(h.dat.F.trace), 'single');
end
h.dat.plot_neu = 0;
h.dat.cl.skewF = skewness(h.dat.F.trace,0,2)';


F=cat(2,h.dat.F.Fcell{:});
zF = zscore(F,0,2);
d = 10;
dzF = zscore(F(:,(1+d):end)-F(:,1:end-d),0,2);
% cellid = find(h.dat.cl.iscell); % to visualize correlation, 
% h.dat.cl.C_of_zF=zF(cellind,:)*zF(cellind,:)'/size(zF,2);
% would be easier to see.  
h.dat.cl.C_of_zF=zF*zF'/size(zF,2);
h.dat.cl.C_of_dzF=dzF*dzF'/size(dzF,2);
h.dat.cl.IDX = recalc_IDX(h);



function [cell_index,roi_index]=get_correlated_roi(C_of_zF,cell_id,threshold)
%  [cell_index,roi_index]=get_correlated_roi(C_of_zF,cell_id,threshold);
%  Input: 
%  C_of_zF:  cross-covariance or normalized cross-correlation of F. The value is between -1 to 1. 
%  cell_id:  cell id within iscell=1 group. 
%  threshold: Correlation higher than this value is returned.

cell_index=find(C_of_dzF(cell_id,:)>threshold);




function pushbutton61_Callback(hObject, eventdata, h)
% keep TOP variance region
h.dat.cl.topregion = zeros(h.dat.cl.Ly, h.dat.cl.Lx);
xs = repmat(1:h.dat.cl.Lx, h.dat.cl.Ly, 1);
ys = repmat((1:h.dat.cl.Ly)', 1, h.dat.cl.Lx);

for k = 1:length(h.dat.stat)
    if ~isempty(h.dat.stat(k).Vregion)
        [~, itop] = max(h.dat.stat(k).Vregion);
        h.dat.cl.topregion(h.dat.stat(k).region{itop}) = 1;
        h.dat.cl.npix(k) = h.dat.stat(k).npixels(itop);
    end
end

h = buildLambdaValue(h);
h=redraw_figure(h);
guidata(hObject,h)

function pushbutton18_Callback(hObject, eventdata, h)
function slider5_Callback(hObject, eventdata, h)
function slider5_CreateFcn(hObject, eventdata, h)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function slider4_Callback(hObject, eventdata, h)
function slider4_CreateFcn(hObject, eventdata, h)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function slider2_Callback(hObject, eventdata, h)
function slider2_CreateFcn(hObject, eventdata, h)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
    
    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end
    
function slider1_CreateFcn(hObject, eventdata, h)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% function pb_binarymask_Callback(hObject, eventdata, h)
% % binary mask
% h.dat.img0.V = ones(h.dat.cl.Ly, h.dat.cl.Lx);
% h = buildLambdaValue(h);
% 
% guidata(hObject,h);
% h=redraw_figure(h);


% function pb_variancemask_Callback(hObject, eventdata, h)
% % variance explained mask
% h.dat.img0.V = reshape(h.dat.res.M, h.dat.cl.Ly, h.dat.cl.Lx)/h.dat.cl.MeanM;
% h = buildLambdaValue(h);
% guidata(hObject,h);
% h=redraw_figure(h);


function pb_defaulthue_Callback(hObject, eventdata, h)
% original default hue
h.dat.cl.rands   = h.dat.cl.rands_orig;
h.dat.img1.H       = reshape(h.dat.cl.rands(h.dat.res.iclust), h.dat.cl.Ly, h.dat.cl.Lx);
h.dat.img2.H       = reshape(h.dat.cl.rands(h.dat.res.iclust), h.dat.cl.Ly, h.dat.cl.Lx);
guidata(hObject,h);
h=redraw_figure(h);

function pb_randomizehue_Callback(hObject, eventdata, h)
% randomize hue
rng('shuffle') 
unique_iclust = unique(h.dat.res.iclust(:));
h.dat.cl.rands     = rand(1, h.dat.ops.Nk);
h.dat.cl.rands(1)  = .15;
h.dat.img1.H       = reshape(h.dat.cl.rands(h.dat.res.iclust), h.dat.cl.Ly, h.dat.cl.Lx);
h.dat.img2.H       = reshape(h.dat.cl.rands(h.dat.res.iclust), h.dat.cl.Ly, h.dat.cl.Lx);

guidata(hObject,h);
h=redraw_figure(h);



function edit33_Callback(hObject, eventdata, h)
h.dat.cl.pixthresh_percent = str2double(get(h.edit33,'String'));
h = exclude_pixels_percent(h);
h = buildLambdaValue(h);
h=redraw_figure(h);
guidata(hObject,h);

function h = exclude_pixels_percent(h)
h.dat.cl.excl_pix_perc = zeros(h.dat.cl.Ly, h.dat.cl.Lx);
for k = 1:h.dat.ops.Nk
    which_pix = find(h.dat.res.iclust==k);
   [Msort, isort] = sort(h.dat.res.M(which_pix), 'ascend'); 
   Msort = cumsum(Msort);
   Msort = 100 * Msort/max(Msort);
   ifi = find(Msort>h.dat.cl.pixthresh_percent, 1);
   h.dat.cl.excl_pix_perc(which_pix(isort(1:ifi))) = 1;
end

function edit33_CreateFcn(hObject, eventdata, h)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit34_Callback(hObject, eventdata, h)
function edit34_CreateFcn(hObject, eventdata, h)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function pushbutton67_Callback(hObject, eventdata, h)
function edit37_Callback(hObject, eventdata, h)
function edit37_CreateFcn(hObject, eventdata, h)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function pushbutton70_Callback(hObject, eventdata, h)
function pushbutton65_Callback(hObject, eventdata, h)


function edit36_Callback(hObject, eventdata, h)
function edit36_CreateFcn(hObject, eventdata, h)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_Compactness_Callback(hObject, eventdata, h)
h.dat.res.Mrs_thresh = str2double(get(h.edit_Compactness,'String'));
h = splitROIleftright(h);
h = buildLambdaValue(h);
h=redraw_figure(h);
guidata(hObject,h);

function edit_Compactness_CreateFcn(hObject, eventdata, h)
set(hObject,'String', num2str(2));
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function pb_saveprocfile_Callback(hObject, eventdata, h)
filename = [h.dat.filename(1:end-4) '_proc.mat'];
fprintf('Saving results %s ...\n',filename)
h.dat.F.trace = [];

dat = h.dat;
if isfield(dat,'reg_data'),dat=rmfield(dat,'reg_data'); end % to remove huge data. 

try
    save(filename, '-struct', 'dat')
catch
    save(filename,'-v7.3', 'dat')
end

printname = [h.dat.filename(1:end-4) '.png'];
fprintf('Printing %s ...\n',printname);
set(gcf,'PaperPositionMode','auto');
print(gcf,'-dpng',printname)


function pushbutton79_Callback(hObject, eventdata, h)
function pushbutton80_Callback(hObject, eventdata, h)

function edit_PixelCountHigh_Callback(hObject, eventdata, h)
h.dat.cl.npix_high = str2double(get(h.edit_PixelCountHigh,'String'));
h = splitROIleftright(h);
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
h = splitROIleftright(h);
h = buildLambdaValue(h);
h=redraw_figure(h);
guidata(hObject,h);

function edit_PixelCountLow_CreateFcn(hObject, eventdata, h)
set(hObject,'String', num2str(30));
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupmenu11_Callback(hObject, eventdata, h)
function popupmenu11_CreateFcn(hObject, eventdata, h)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit41_Callback(hObject, eventdata, h)
function edit41_CreateFcn(hObject, eventdata, h)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit42_Callback(hObject, eventdata, h)
function edit42_CreateFcn(hObject, eventdata, h)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function figure1_ResizeFcn(hObject, eventdata, h)

function edit43_Callback(hObject, eventdata, h)
h.dat.cl.pixthresh_var = str2double(get(h.edit43,'String'));
h = excluded_pixels(h);
h = buildLambdaValue(h);
h=redraw_figure(h);
guidata(hObject,h);


function edit43_CreateFcn(hObject, eventdata, h)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit44_Callback(hObject, eventdata, h)
function edit44_CreateFcn(hObject, eventdata, h)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

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
        if h.quadvalue(j,i)==1
            set(h.(sprintf('Q%d%d', j,i)), 'BackgroundColor','yellow'); 
        end
    end
end
set(h.(sprintf('Q%d%d', iy,ix)), 'BackgroundColor','red'); 

% --- Executes on button press in full.
function full_Callback(hObject, eventdata, h)
h.dat.ylim = [0 h.dat.cl.Ly];
h.dat.xlim = [0 h.dat.cl.Lx];
guidata(hObject,h);
h=redraw_figure(h);

function quadrant(hObject, h, iy, ix)
h.dat.ylim = [h.dat.figure.y0all(iy) h.dat.figure.y1all(iy+1)];
h.dat.xlim = [h.dat.figure.x0all(ix) h.dat.figure.x1all(ix+1)];
h.quadvalue(iy, ix) = 1;

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
    case 'f'
%         % flip currently selected unit
%         h.dat.cl.manual(h.dat.F.ichosen) = .5 - h.dat.cl.iscell(h.dat.F.ichosen);
%         h = splitROIleftright(h);
%         h = buildLambdaValue(h);
%         guidata(hObject,h);
%         if h.dat.maxmap==1
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
        h.dat.map = 1;
        h=pb_rois_display_Callback(hObject, eventdata, h);
        last_keypress = eventdata.Key;
    case 'c'
         set(h.uipanel_Left,'Title','Correlated ROIs');
        h.dat.map = 2;
        h=pb_CorrROI_display_Callback(hObject, eventdata, h);
        last_keypress = eventdata.Key;
        
    case 'w'
        set(h.uipanel_Left,'Title','weighted image of corr-pixel ("p" to proc on/off)');
        h.dat.map = 2;
        h=pb_meandisplay_Callback(hObject, eventdata, h);
        last_keypress = eventdata.Key;
  
    case 'e'
       
        h.dat.map = 3;
        if h.dat.maxmap>2
            pb_Rede_display_Callback(hObject, eventdata, h);
        end
        last_keypress = eventdata.Key;
    case 'r'
        h.dat.map = 3;
        if h.dat.maxmap>3
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
x = min(max(1, round(x)), h.dat.cl.Lx);
y = min(max(1, round(y)), h.dat.cl.Ly);

h.dat.F.ichosen = h.dat.res.iclust(y, x);
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
%             h.dat.cl.manual(h.dat.F.ichosen) = .5 - h.dat.cl.iscell(h.dat.F.ichosen);
%             h = splitROIleftright(h);
%             h = buildLambdaValue(h);
        case 'open' % double click
            % unpin the manual selection on this cell
%             h.dat.cl.manual(h.dat.F.ichosen) = 0;
%             h = splitROIleftright(h);
%             h = buildLambdaValue(h);
        case 'extend' % shift left
%             h.dat.cl.redcell(h.dat.F.ichosen) = 1 -  h.dat.cl.redcell(h.dat.F.ichosen);
%             
%             if h.dat.cl.redcell(h.dat.F.ichosen)==1
%                 h.dat.cl.rands(h.dat.F.ichosen) = 0;
%             else
%                 h.dat.cl.rands(h.dat.F.ichosen) = h.dat.cl.rands_orig(h.dat.F.ichosen);
%             end
%             h.dat.img1.H       = reshape(h.dat.cl.rands(h.dat.res.iclust), h.dat.cl.Ly, h.dat.cl.Lx);
%             h.dat.img2.H       = reshape(h.dat.cl.rands(h.dat.res.iclust), h.dat.cl.Ly, h.dat.cl.Lx);
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
            CallbackFcn=@(hObject,eventdata)RoiGui_004('enter_svm_value', hObject,eventdata,h);            
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
    h.dat.cl.iscell(h.dat.F.ichosen),x,y);
fprintf('Compactness=%2.1f, npix=%d\n', h.dat.cl.Mrs(h.dat.F.ichosen),...
    h.dat.cl.npix(h.dat.F.ichosen));
fprintf('Skew(F)=%3.3f, ',h.dat.cl.skewF(h.dat.F.ichosen));
fprintf('Eccentricity=%3.3f, ',h.dat.cl.Eccentricity(h.dat.F.ichosen));
fprintf('Solidity=%3.3f, ',h.dat.cl.Solidity(h.dat.F.ichosen));
fprintf('meanV=%3.3f,',h.dat.cl.V(h.dat.F.ichosen)/h.dat.cl.npix(h.dat.F.ichosen));
fprintf('manual label=%3.3d\n',h.dat.cl.manual_type(h.dat.F.ichosen));

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
x = min(max(1, round(x)), h.dat.cl.Lx);
y = min(max(1, round(y)), h.dat.cl.Ly);


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
x1 = min(max(1, round(x1)), h.dat.cl.Lx);
y1 = min(max(1, round(y1)), h.dat.cl.Ly);
hold on; plotH=plot(x1,y1,'ro');

% another center of k-means clustering
x2 = ceil(ginput(1));
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
    fprintf('Loading svd file: %s\n', SVDFile);
    
    if exist(SVDFile,'file')    
        waitH = waitbar(0.2,sprintf('Loding %s',SVDFile));
        h.dat.svd=   load(SVDFile, 'U','Sv');
        waitbar(1,waitH);
        delete(waitH);
    else
        error('SVD file %s NOT FOUND.',SVDFile);
    end
    
end

L=h.dat.res.iclust;

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

% figH=myfigure('Split or not');
% set(figH,'Position',[400  100  360  640]);
% clf;subplot(2,1,1);
Sat=reshape(h.dat.res.probabilities,h.dat.res.Ly,h.dat.res.Lx);
% IMG=hsv2rgb(cat(3,h.dat.img1.H,Sat,Val));
% imagesc(IMG(J,I,:));
% drawnow;

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
%  subplot(2,1,2);cla;
%  for ii=1:Nlabel
%      pind=find(idx==ii);
%      plotH(ii)=plot3(score(pind,1),score(pind,2),score(pind,3),[Col{ii},'.']);hold on;
%  end
% grid on;

% recover the divided image 

% subplot(2,1,1);cla;

Col = [0.9,0.4,0.6];
% Val = L~=0;

Hue = h.dat.img1.H;
for ii=1:Nlabel
    original_ind = ind(idx==ii);
    Hue(original_ind)=Col(ii);
%     Sat(original_ind)=1;
end

IMG=hsv2rgb(cat(3,Hue,Sat,Val));
% image(IMG(J,I,:));

SplitOrNot = ROISplit_GUI('title','ROI split',...
    'CData',IMG(J,I,:),...
    'PlotData',score,...
    'PlotData_Class',idx)

unique_labels = unique(h.dat.res.iclust1);
labelN = max(unique_labels);
switch SplitOrNot
    case 'Split'
        for ii=1:Nlabel
            BW=zeros(size(h.dat.res.iclust));
            BW(ind(idx==ii))=1;
            CC = bwconncomp(BW);
            [~,max_id] = max(cellfun(@length,CC.PixelIdxList));
            for jj=1:length(CC.PixelIdxList)
                 h.dat.res.iclust1(CC.PixelIdxList{jj})=labelN+1;
                 labelN= labelN+1;
            end          
           
        end
    case 'No, I change my mind'
        
    otherwise
        error('Unknown answer %s',SplitOrNot)
end
guidata(hObject,h);

 % cell selection function
function CellSelectionButtonDownFcn(hObject,eventdata,h)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pos=get(gca,'CurrentPoint');

x = round(pos(1,1));
y  = round(pos(1,2));
x = min(max(1, round(x)), h.dat.cl.Lx);
y = min(max(1, round(y)), h.dat.cl.Ly);

if isempty(h.dat.F.ichosen_append), h.dat.F.ichosen_append=[];end
h.dat.F.ichosen = h.dat.res.iclust(y, x);
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
            h.dat.cl.manual(h.dat.F.ichosen) = .5 - h.dat.cl.iscell(h.dat.F.ichosen);
            h = splitROIleftright(h);
            h = buildLambdaValue(h);
        case 'open' % double click
            % unpin the manual selection on this cell
            h.dat.cl.manual(h.dat.F.ichosen) = 0;
            h = splitROIleftright(h);
            h = buildLambdaValue(h);
        case 'extend' % shift left
            h.dat.cl.redcell(h.dat.F.ichosen) = 1 -  h.dat.cl.redcell(h.dat.F.ichosen);
            
            if h.dat.cl.redcell(h.dat.F.ichosen)==1
                h.dat.cl.rands(h.dat.F.ichosen) = 0;
            else
                h.dat.cl.rands(h.dat.F.ichosen) = h.dat.cl.rands_orig(h.dat.F.ichosen);
            end
            h.dat.img1.H       = reshape(h.dat.cl.rands(h.dat.res.iclust), h.dat.cl.Ly, h.dat.cl.Lx);
            h.dat.img2.H       = reshape(h.dat.cl.rands(h.dat.res.iclust), h.dat.cl.Ly, h.dat.cl.Lx);
            
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
    h.dat.cl.manual(h.dat.F.ichosen),...
    h.dat.cl.iscell(h.dat.F.ichosen),x,y);
fprintf('Compactness=%2.1f, npix=%d\n', h.dat.cl.Mrs(h.dat.F.ichosen),...
    h.dat.cl.npix(h.dat.F.ichosen));
fprintf('Skew(F)=%3.3f, ',h.dat.cl.skewF(h.dat.F.ichosen));
fprintf('Eccentricity=%3.3f, ',h.dat.cl.Eccentricity(h.dat.F.ichosen));
fprintf('Solidity=%3.3f, ',h.dat.cl.Solidity(h.dat.F.ichosen));
fprintf('meanV=%3.3f,',h.dat.cl.V(h.dat.F.ichosen)/h.dat.cl.npix(h.dat.F.ichosen));
fprintf('manual label=%d\n',h.dat.cl.manual_type(h.dat.F.ichosen));


h=update_figure(h);
guidata(hObject,h);

function h=update_figure(h)
Mode = h.ModeSelectionButtonGroup.SelectedObject.Tag;

switch Mode
    case {'tg_DisplayMode','tg_ROISplitMode'}
        PlotHue = h.dat.cl.rands(h.dat.F.ichosen_append);
        PlotColor =permute(cat(2,PlotHue(:),ones(length(PlotHue),2)),[1,3,2]);
        PlotColor = hsv2rgb(PlotColor);
        h=redraw_fluorescence_multi(h,PlotColor);
        
        h = buildSat(h);
        h = buildLambdaValue(h);
        h = buildHue(h);
        
        %fprintf('calculating cell mask\n')
%         h.dat.cellmask=zeros(h.dat.cl.Ly, h.dat.cl.Lx);
%         for ii=1:length(h.dat.F.ichosen_append)
% %             h.dat.cellmask(h.dat.res.iclust==h.dat.F.ichosen_append(ii))=1;
%             h.dat.cellmask(h.dat.stat(h.dat.F.ichosen_append(ii)).ipix_edge)=1;
%         end
%         h.dat.cellmask = edge(h.dat.cellmask);
        h=redraw_figure(h);
        
    case 'tg_EllipseMode'
%         h=add_imellipse_Callback(hObject, eventdata, h);
    case 'tg_ROIMergeMode'
        h=redraw_corrimg(h);
    case 'tb_SVMMode'
        PlotHue = h.dat.cl.rands(h.dat.F.ichosen_append);
        PlotColor =permute(cat(2,PlotHue(:),ones(length(PlotHue),2)),[1,3,2]);
        PlotColor = hsv2rgb(PlotColor);
        h=redraw_fluorescence_multi(h,PlotColor);
        
        h = buildHue(h);
        h = buildSat(h);
        h = buildLambdaValue(h);
    
        
        %fprintf('calculating cell mask\n')
%         h.dat.cellmask=zeros(h.dat.cl.Ly, h.dat.cl.Lx);
%         for ii=1:length(h.dat.F.ichosen_append)
%             h.dat.cellmask(h.dat.res.iclust==h.dat.F.ichosen_append(ii))=1;
%         end
%         h.dat.cellmask = edge(h.dat.cellmask);
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
x = min(max(1, round(x)), h.dat.cl.Lx);
y = min(max(1, round(y)), h.dat.cl.Ly);

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

if h.dat.map>1
    h=redraw_meanimg(h);
end

% if h.dat.map < h.dat.maxmap
%     h.dat.map = h.dat.map + 1;
%     h=redraw_meanimg(h);
% else
%     h.dat.map = 1;
%     h=redraw_figure(h);
% end
guidata(hObject,h);


% --- Executes on button press in pb_rois_display.
function h=pb_rois_display_Callback(hObject, eventdata, h)
h.dat.map = 1;
h=redraw_figure(h);
guidata(hObject,h);


% --- Executes on button press in pb_meandisplay.
function h=pb_meandisplay_Callback(hObject, eventdata, h)
h.dat.map = 2;
h=redraw_meanimg(h);
guidata(hObject,h);

% --- Executes on button press in pb_Rede_display.
function pb_Rede_display_Callback(hObject, eventdata, h)
 h.dat.map = 3;
h=redraw_meanimg(h);
guidata(hObject,h);

% RED CORRECTED BUTTON
function pb_redcorr_display_Callback(hObject, eventdata, h)
h.dat.map = 4;
h=redraw_meanimg(h);
guidata(hObject,h);






% --- Executes on button press in tg_showneuropil.
function tg_showneuropil_Callback(hObject, eventdata, h)
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
if isempty(h.dat.F.ichosen_append)
    x_trace = [1:NT];
    y_trace=my_conv_local(double(h.dat.F.trace(h.dat.F.ichosen,:)), 3);
else
    x_trace = [1:NT]; %repmat([1:NT NaN],size(h.dat.F.ichosen_append,1),1)';
    y_trace=zscore(my_conv_local(double(h.dat.F.trace(h.dat.F.ichosen_append,:)), 3),0,2);
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
fprintf('Skew(F)=%3.3f\n',h.dat.cl.skewF(h.dat.F.ichosen));
y_neurop=my_conv_local(medfilt1(double(h.dat.F.neurop(h.dat.F.ichosen,:)), 3), 3);

if get(h.tg_showneuropil,'Value')
    if isfield(h.dat.F, 'neurop')
      
        set(h.plot_fluorescence_baseF,'XData',1:NT,'YData',y_neurop,'Color','c','LineStyle','--');
    else
        set(h.plot_fluorescence_baseF,'XData',NaN,'YData',NaN);

    end
else
    set(h.plot_fluorescence_baseF,'XData',NaN,'YData',NaN);
end
% y = double(h.dat.F.trace(h.dat.F.ichosen,:))';
x = [ones(length(y_neurop),1), y_neurop'];
% b=regress(y_trace,x);
b = x'*x\x'*y_trace';
newydata=y_trace-b(2)*y_neurop-b(1);
newydata=my_conv_local(newydata',3)+mean(y_trace);

set(h.plot_fluorescence_FCell,'XData',1:NT,'YData',newydata,'Color','r');
max_y = max([max_y;newydata]);
min_y = min([min_y;newydata]);
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
    I = hsv2rgb(cat(3, h.dat.img1.H, h.dat.img1.Sat, h.dat.img1.V));
else % show UN-selected one
    I = hsv2rgb(cat(3, h.dat.img2.H, h.dat.img2.Sat, h.dat.img2.V));
end

I = min(I, 1);
% axes(h.axes_left); imagesc(I);
set(h.left_imageH,'CData',I);
set(h.axes_left,'XLim',[h.dat.xlim],'YLim',[h.dat.ylim]);



function h=redraw_meanimg(h)

if h.dat.procmap
    I = h.dat.mimg_proc(:,:,h.dat.map);
else    
    I = h.dat.mimg(:,:,h.dat.map);
end

mu = median(I(:));
sd1 = mean(abs(I(I<mu) - mu));
sd2 = mean(abs(I(I>mu) - mu));

% axes(h.axes_left); imagesc(I, mu + 5*[-sd1 sd2]);
set(h.left_imageH,'CData',I);
% xlim([h.dat.xlim]); ylim([h.dat.ylim]);
set(h.axes_left,'XLim',[h.dat.xlim],'YLim',[h.dat.ylim],'CLim',mu + [-2.5*sd1 4*sd2]);

% axes(h.axes_right); imagesc(I, mu + 5*[-sd1 sd2]);
% xlim([h.dat.xlim]); ylim([h.dat.ylim]);
% set(h.left_imageH,'CData',I);
% xlim([h.dat.xlim]); ylim([h.dat.ylim]);
% set(h.axes_right,'XLim',[h.dat.xlim],'YLim',[h.dat.ylim],'CLim',mu + [-2.5*sd1 4*sd2]);


drawnow

function h=show_corr_roi(h)
%%
IDX=h.dat.cl.IDX;
unique_cluster = unique(IDX);
unique_cluster(unique_cluster==0)=[];



% roi_samecell=find(h.dat.cl.C_of_zF(h.dat.F.ichosen,:)>threshold);


Vmap = h.dat.img1.V; % value, 1 is brighter


figh=myfigure('CorrROI');clf;% set(figh,'Position',);
for ii=unique_cluster(:)'
    subplot(2,1,1);
    roi_samecell = find(IDX==ii);
    Hmap = zeros(1,length(h.dat.res.iclust)); % hue color map, make the same cell same color.
    Hmap(roi_samecell)=1;
    Hmap       = reshape(Hmap(h.dat.res.iclust), h.dat.cl.Ly, h.dat.cl.Lx);
    Smap = zeros(size(h.dat.img1.Sat)); % saturation, 1 is more pure color, 0 is white.
    for jj=roi_samecell(:)'
        Smap(h.dat.res.iclust==jj)=1;
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

selected_region = edge(h.dat.res.iclust==h.dat.F.ichosen);

Smap = zeros(size(h.dat.img1.Sat)); % saturation, 1 is more pure color, 0 is white.
% color the selected cell, and let the others remain white. 
for ii=1:length(roi_samecell)
    Smap(h.dat.res.iclust==roi_samecell(ii))=1;
end

Hmap      = reshape(h.dat.cl.rands(h.dat.res.iclust), h.dat.cl.Ly, h.dat.cl.Lx);

% Hmap = zeros(1,length(h.dat.res.iclust)); % hue color map, make the same cell same color.
% Hmap(roi_samecell)=0.6; % blue
% Hmap=h.dat.cl.rands;
% Hmap       = reshape(Hmap(h.dat.res.iclust), h.dat.cl.Ly, h.dat.cl.Lx);


Vmap = h.dat.img1.V+selected_region; % value, 1 is brighter



I = hsv2rgb(cat(3, Hmap, Smap, Vmap));
I = min(I, 1);
% axes(h.axes_left); imagesc(I);
set(h.left_imageH,'CData',I);
set(h.axes_left,'XLim',[h.dat.xlim],'YLim',[h.dat.ylim]);
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
            h=pb_loadmovie_Callback(hObject, eventdata, h);
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

function h = splitROIleftright(h)

contents=get(h.popup_FilterSelection,'String');
SelectedFilter = contents{get(h.popup_FilterSelection,'Value')}


switch SelectedFilter
    case 'Param filter'
        h.dat.cl.iscell = h.dat.cl.Mrs <h.dat.res.Mrs_thresh & ...
            h.dat.cl.npix <h.dat.cl.npix_high & ...
            h.dat.cl.npix >h.dat.cl.npix_low;
        
        % h.dat.cl.iscell = h.dat.cl.iscell  & (h.dat.cl.nreg      <h.dat.cl.nreg_max);
        % h.dat.cl.iscell = h.dat.cl.iscell  & (h.dat.cl.npix_par  <h.dat.cl.npix_par_max);
        % h.dat.cl.iscell = h.dat.cl.iscell  & (h.dat.cl.npix_res  <h.dat.cl.npix_res_max);
        % h.dat.cl.iscell = h.dat.cl.iscell  & (h.dat.cl.mrs_parent<h.dat.cl.mrs_parent_max);
        % h.dat.cl.iscell = h.dat.cl.iscell  & (h.dat.cl.VperPix   >h.dat.cl.VperPix_min);
        h.dat.cl.skewF_low = str2double(get(h.edit_minSkewF,'String'));
        h.dat.cl.EccentricityMax = str2double(get(h.edit_EccentricityMax,'String'));
        h.dat.cl.SolidityMin = str2double(get(h.edit_SolidityMin,'String'));
        
        h.dat.cl.iscell = h.dat.cl.iscell  & (h.dat.cl.skewF(:)'     >h.dat.cl.skewF_low);
        h.dat.cl.iscell = h.dat.cl.iscell  & (h.dat.cl.Eccentricity(:)'     <h.dat.cl.EccentricityMax);
        h.dat.cl.iscell = h.dat.cl.iscell  & (h.dat.cl.Solidity(:)'     >h.dat.cl.SolidityMin);
        
       
    case 'SVM filter'
        h.dat.cl.iscell = zeros(size(h.dat.cl.Mrs));
        
        if (h.cb_SVM_CellBody.Value)
            h.dat.cl.iscell(h.dat.cl.predicted_type==1)=1;
        end
        
        if (h.cb_SVM_Dendrite.Value)
            h.dat.cl.iscell(h.dat.cl.predicted_type==2)=1;
        end
        
        if (h.cb_SVM_Axon.Value)
            h.dat.cl.iscell(h.dat.cl.predicted_type==3)=1;
        end
        
         h.dat.cl.iscell(h.dat.cl.predicted_type==9)=0;
        
        
    case 'none'
        
    otherwise
        error('Unknown filter %s',SelectedFilter);
end

h.dat.cl.iscell(max(h.dat.res.iclust(:)))=0;
h.dat.cl.iscell = double(h.dat.cl.iscell);

% always manual selections can overwrite. 
h.dat.cl.iscell(h.dat.cl.manual>1e-3) = 1;
h.dat.cl.iscell(h.dat.cl.manual<-1e-3) = 0;
h.dat.cl.k1 = reshape(h.dat.cl.iscell(h.dat.res.iclust), h.dat.cl.Ly, h.dat.cl.Lx);

function h = setOriginalThresh(h)
h.dat.res.Mrs_thresh = h.dat.res.Mrs_thresh_orig;
h.dat.cl.npix_low = h.dat.cl.npix_low_orig;
h.dat.cl.npix_high = h.dat.cl.npix_high_orig;

function h = buildHue(h)
Tag=h.ModeSelectionButtonGroup.SelectedObject.Tag;
h.dat.img1.H       = reshape(h.dat.cl.rands(h.dat.res.iclust), h.dat.cl.Ly, h.dat.cl.Lx);
h.dat.img2.H       = h.dat.img1.H;
        
switch Tag
    case {'tg_DisplayMode','tg_EllipseMode','tg_ROIMergeMode','tg_ROISplitMode'}    
     % do nothing
    case 'tb_SVMMode'
        
       for ii=1:h.dat.ops.Nk
            %             h.dat.cellmask(h.dat.res.iclust==h.dat.F.ichosen_append(ii))=1;
            edge_pix = h.dat.stat(ii).ipix_edge;
            type_index = find([h.dat.cl.type_number_table{:,1}]==h.dat.cl.manual_type(ii));
            HSV=rgb2hsv(h.dat.cl.type_number_table{type_index,3}); 
            if type_index~=1
                h.dat.img1.H(edge_pix)=HSV(1);
                h.dat.img2.H(edge_pix)=HSV(1);
            end
        end
    otherwise
        error('Unknown Tag %s',Tag);
end


function h = buildSat(h) % control gray to color mode. 0 is gray. 1 is color.

Tag=h.ModeSelectionButtonGroup.SelectedObject.Tag;

switch Tag
    case {'tg_DisplayMode','tg_EllipseMode','tg_ROIMergeMode','tg_ROISplitMode'}
        Sat = ones(h.dat.cl.Ly, h.dat.cl.Lx);
        % make background gray
        Sat(h.dat.res.iclust==max(h.dat.res.iclust(:)))=0;
        
        % make selected cell white fringed.
        h.dat.cellmask=zeros(h.dat.cl.Ly, h.dat.cl.Lx);
        for ii=1:length(h.dat.F.ichosen_append)
%             h.dat.cellmask(h.dat.res.iclust==h.dat.F.ichosen_append(ii))=1;
            edge_pix = h.dat.stat(h.dat.F.ichosen_append(ii)).ipix_edge;
            h.dat.cellmask(edge_pix)=1;
            Sat(edge_pix)=0;
        end
           
        h.dat.img1.Sat     = Sat;
        h.dat.img2.Sat     = Sat;
   
    case 'tb_SVMMode'
        
        Sat = ones(h.dat.cl.Ly, h.dat.cl.Lx);
        % make background faintly colored 
        Sat(h.dat.res.iclust==max(h.dat.res.iclust(:)))=0;
        
        for ii=1:h.dat.ops.Nk
            %             h.dat.cellmask(h.dat.res.iclust==h.dat.F.ichosen_append(ii))=1;
            edge_pix = h.dat.stat(ii).ipix_edge;
            type_index = find([h.dat.cl.type_number_table{:,1}]==h.dat.cl.manual_type(ii));
            HSV=rgb2hsv(h.dat.cl.type_number_table{type_index,3}); 
            if type_index~=1
                Sat(edge_pix)=HSV(2);
            end
        end
        
        h.dat.cellmask=zeros(h.dat.cl.Ly, h.dat.cl.Lx);
        for ii=1:length(h.dat.F.ichosen_append)
            %             h.dat.cellmask(h.dat.res.iclust==h.dat.F.ichosen_append(ii))=1;
            edge_pix = h.dat.stat(h.dat.F.ichosen_append(ii)).ipix_edge;
            h.dat.cellmask(edge_pix)=1;
            Sat(edge_pix)=0;
        end
        
        h.dat.img1.Sat     = Sat;
        h.dat.img2.Sat     = Sat;
        
    otherwise
        error('Unkown mode %s',Tag);
end


function h = buildLambdaValue(h) % V of HSV color; it controls Brightness. 
% h.dat.img1.V       = h.dat.img0.V .* h.dat.cl.k1 .* (1-h.dat.cl.excluded_pixels)...
%     .*(1-h.dat.cl.excluded_regions).*(1-h.dat.cl.excl_pix_perc).*h.dat.cl.topregion;
% h.dat.img2.V       = h.dat.img0.V .* (~h.dat.cl.k1);
Tag=h.ModeSelectionButtonGroup.SelectedObject.Tag;

BGval = get(h.slider_BGLevel,'Value');
BG=BGval*h.dat.img0.BackGround.*h.dat.res.M0;

Brightness = get(h.slider_Sat,'Value');
h.dat.img1.V =BG+Brightness^2*h.dat.res.M0.*~h.dat.img0.BackGround...
    .* h.dat.cl.k1 .* (1-h.dat.cl.excluded_pixels)...
    .*(1-h.dat.cl.excluded_regions)...
    .*(1-h.dat.cl.excl_pix_perc).*h.dat.cl.topregion;
h.dat.img2.V =BG+Brightness^2*h.dat.res.M0.*(~h.dat.cl.k1.*~h.dat.img0.BackGround);
       
switch Tag
   case {'tg_DisplayMode','tg_EllipseMode','tg_ROIMergeMode','tg_ROISplitMode'}
    
      % do nothing more. 
               
    case 'tb_SVMMode'
        
        for ii=1:h.dat.ops.Nk
            %             h.dat.cellmask(h.dat.res.iclust==h.dat.F.ichosen_append(ii))=1;
            edge_pix = h.dat.stat(ii).ipix_edge;
            type_index = find([h.dat.cl.type_number_table{:,1}]==h.dat.cl.manual_type(ii));
            HSV=rgb2hsv(h.dat.cl.type_number_table{type_index,3});
            
            if type_index ~=1
                  h.dat.img1.V(edge_pix)=HSV(3);
            end
          
%             h.dat.img2.V(edge_pix)=HSV(3);
        end
        
    otherwise
        error('Unknown Tag %s',Tag);
end
    
% iselect = h.dat.res.iclust==h.dat.F.ichosen;
% h.dat.img1.V(iselect) = h.dat.cl.k1(iselect);
% h.dat.img2.V(iselect) = ~h.dat.cl.k1(iselect);


for ii=h.dat.F.ichosen_append(:)'
    edge_pix = h.dat.stat(ii).ipix_edge;
   
    h.dat.img1.V(edge_pix)=1;
    h.dat.img2.V(edge_pix)=1;
end

function h = excluded_pixels(h)
h.dat.cl.excluded_pixels = h.dat.res.M<h.dat.cl.pixthresh_var;
h.dat.cl.excluded_pixels = reshape(h.dat.cl.excluded_pixels, h.dat.cl.Ly, h.dat.cl.Lx);


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
function pb_movie_display_Callback(hObject, eventdata, handles)
% hObject    handle to pb_movie_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h.dat.map = 5;
h=redraw_meanimg(h);
guidata(hObject,h);




% --- Executes on button press in pb_loadmovie.
function handles=pb_loadmovie_Callback(hObject, eventdata, handles)
% hObject    handle to pb_loadmovie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

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

if FileExists(fullfile(RegFastFilePath,RegFastFile))
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
ind_i = 1:max(handles.dat.ylim);
ind_j = 1:max(handles.dat.xlim);

handles.dat.reg_fast_filename=fullfile(RegFastFilePath,RegFastFile);

% from RoiGui_004, use mlttiff class to access multiple Tiff files as a single object. 
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
h = splitROIleftright(h);
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

function merge_two_roi(hObject,eventdata,h)


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
     handles=pb_loadmovie_Callback(hObject, eventdata, handles);
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

h.dat.map = 1;
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


% --- Executes on button press in tg_ROIMergeMode.
function tg_ROIMergeMode_Callback(hObject, eventdata, handles)
% hObject    handle to tg_ROIMergeMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tg_ROIMergeMode


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
            
        ButtonDownFcn=@(hObject,eventdata)RoiGui_004('CellSelectionButtonDownFcn',...
            hObject,eventdata,guidata(hObject));
        set(handles.left_imageH,'ButtonDownFcn',ButtonDownFcn);
         set(handles.uipanel_Left,'Title','DisplayMode');
    case 'tg_EllipseMode'
        
        ButtonDownFcn=@(hObject,eventdata)RoiGui_004('EllipseSelectionButtonDownFcn',...
            hObject,eventdata,guidata(hObject));
        set(handles.left_imageH,'ButtonDownFcn',ButtonDownFcn);
         set(handles.uipanel_Left,'Title','EllipseMode');
    case 'tg_ROIMergeMode'
        
        ButtonDownFcn=@(hObject,eventdata)RoiGui_004('CellSelectionButtonDownFcn',...
            hObject,eventdata,guidata(hObject));
        set(handles.left_imageH,'ButtonDownFcn',ButtonDownFcn);
        
        handles=renew_cluster(handles);
    case 'tg_ROISplitMode'
        
        ButtonDownFcn=@(hObject,eventdata)RoiGui_004('CellSplitButtonDownFcn',...
            hObject,eventdata,guidata(hObject));
        set(handles.left_imageH,'ButtonDownFcn',ButtonDownFcn);
        set(handles.uipanel_Left,'Title','ROISplitMode');
    case 'tb_SVMMode'
        
        ButtonDownFcn=@(hObject,eventdata)RoiGui_004('CellSelectionButtonDownFcn',...
            hObject,eventdata,guidata(hObject));
%         ButtonDownFcn=@(hObject,eventdata)RoiGui_004('SVMSelectionButtonDownFcn',...
%             hObject,eventdata,guidata(hObject));
        set(handles.left_imageH,'ButtonDownFcn',ButtonDownFcn);
         set(handles.uipanel_Left,'Title','SVMMode');
           
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
    [npix_max,npix_max_ind]=  max(h.dat.cl.npix(ind));
    IDXtmp(ind)= ind(npix_max_ind);
end
IDX = IDXtmp;

function h=update_iclust1(h)

IDX=h.dat.cl.IDX;
h.dat.res.iclust1 = h.dat.res.iclust;

% first, assign index to each cluster
h.dat.cl.automerge =IDX;
ind = find(h.dat.cl.automerge==0);
h.dat.cl.automerge(ind)=ind; 

% overwrite the cluster
ind=find(h.dat.cl.manualmerge);
h.dat.cl.roi_cluster = h.dat.cl.automerge;
h.dat.cl.roi_cluster(ind)=h.dat.cl.manualmerge(ind);

unique_cluster=unique(h.dat.cl.roi_cluster);
unique_cluster(unique_cluster==0)=[];

for ii=1:length(unique_cluster)
    ind = find(IDX==unique_cluster(ii));
    [npix_max,npix_max_ind]=  max(h.dat.cl.npix(ind));
    for jj=1:length(ind)
    h.dat.res.iclust1(h.dat.res.iclust1==ind(jj))=ind(npix_max_ind); % set the ID to the ROI with biggest npix
    end
end

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
% Smap = ones(size(h.dat.img1.Sat)); % saturation, 1 is more pure color, 0 is white.
% % Hmap = zeros(1,length(h.dat.res.iclust)); % hue color map, make the same cell same color.
% % Hmap(roi_samecell)=h.dat.cl.rands(h.dat.F.ichosen);
% Hmap       = h.dat.img2.H; reshape(Hmap(h.dat.res.iclust), h.dat.cl.Ly, h.dat.cl.Lx);
% 
% Vmap = h.dat.img2.V; % value, 1 is brighter
% 
% 
% I = hsv2rgb(cat(3, Hmap, Smap, Vmap));
% I = min(I, 1);
% image(I);
Vmap = h.dat.img1.V;
figh=myfigure('CorrROI');clf;% set(figh,'Position',);
for ii=unique_cluster
    subplot(2,1,1);
    roi_samecell = find(h.dat.cl.IDX==ii);
    Hmap = zeros(1,length(h.dat.res.iclust)); % hue color map, make the same cell same color.
    Hmap(roi_samecell)=h.dat.cl.rands(roi_samecell(1)   );
    Hmap       = reshape(Hmap(h.dat.res.iclust), h.dat.cl.Ly, h.dat.cl.Lx);
    Smap = zeros(size(h.dat.img1.Sat)); % saturation, 1 is more pure color, 0 is white.
    for jj=roi_samecell(:)'
        Smap(h.dat.res.iclust==jj)=1;
    end
    I = hsv2rgb(cat(3, Hmap, Smap, Vmap));
    I = min(I,1);
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
        unique_ROI_elements = unique(h.dat.res.iclust);     
        s=cell(length(unique_ROI_elements),1);
        for ii=1:length(unique_ROI_elements)
            if h.dat.cl.iscell(ii)
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
handles = buildLambdaValue(handles) %
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

handles = buildLambdaValue(handles) %
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


% --- Executes on slider movement.
function slider13_Callback(hObject, eventdata, handles)
% hObject    handle to slider13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider13 (see GCBO)
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
h = splitROIleftright(h);
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


h.dat.cl.SolidityMin = str2double(get(h.edit_SolidityMin,'String'));
h = splitROIleftright(h);
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


function h=update_manual_type(hObject,eventdata,h,varargin)

if nargin>=4
    roi_type = varargin{1};
else
    error();
end
% roi_type
% h.dat.F.ichosen_append

h.dat.cl.manual_type(h.dat.F.ichosen_append)=roi_type;
guidata(hObject,h);

% --- Executes on button press in pb_SVM_1_CellBody.
function pb_SVM_1_CellBody_Callback(hObject, eventdata, handles)
% hObject    handle to pb_SVM_1_CellBody (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

update_manual_type(hObject,eventdata,handles,1);

% --- Executes on button press in pb_2_Dendrite.
function pb_2_Dendrite_Callback(hObject, eventdata, handles)
% hObject    handle to pb_2_Dendrite (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
update_manual_type(hObject,eventdata,handles,2);

% --- Executes on button press in pb_3_Axon.
function pb_3_Axon_Callback(hObject, eventdata, handles)
% hObject    handle to pb_3_Axon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
update_manual_type(hObject,eventdata,handles,3);

% --- Executes on button press in pb_SVM_4_Spine.
function pb_SVM_4_Spine_Callback(hObject, eventdata, handles)
% hObject    handle to pb_SVM_4_Spine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
update_manual_type(hObject,eventdata,handles,4);

% --- Executes on button press in pb_SVM_0_Unknown.
function pb_SVM_0_Unknown_Callback(hObject, eventdata, handles)
% hObject    handle to pb_SVM_0_Unknown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
update_manual_type(hObject,eventdata,handles,0);


% --- Executes on button press in pb_SVM_9_Noise.
function pb_SVM_9_Noise_Callback(hObject, eventdata, handles)
% hObject    handle to pb_SVM_9_Noise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

update_manual_type(hObject,eventdata,handles,9);

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

labels = h.dat.cl.manual_type;
teacher_ind = (labels~=0);
% labels = labels(teacher_ind);
unique_classes = unique(labels); unique_classes(unique_classes==0)=[];

data = [h.dat.cl.Mrs; ...
    h.dat.cl.npix;...
    h.dat.cl.Eccentricity;,...
    h.dat.cl.Solidity;...
    h.dat.cl.V;...
    h.dat.cl.skewF];
data = double(data');

% data= double(data(:,teacher_ind));

fprintf('training...\n')

for cc=unique_classes(:)'
    train_label= double(labels==cc);
    
    fprintf('class=%d:(N=%d)\n',cc,nnz(train_label));
    h.dat.svmmodel(cc).param = svmtrain(2*train_label(teacher_ind)-1,...
        data(teacher_ind,:), '-c 400 -g 0.01');
    
    [predict_label, accuracy, dec_values] = svmpredict(randn(h.dat.ops.Nk,1), ...
        data,  h.dat.svmmodel(cc).param); % test the training data
    
    h.dat.cl.predicted_type(predict_label>0)=cc;
end
h = splitROIleftright(h);
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
