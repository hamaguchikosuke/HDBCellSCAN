%% align ROIs from different days 
% 1. load pipeline-processed proc data
%    load the mean data 
%  <To Do>: change this to proc-behav data.

% 
% file{1}='E:\home\ImagingData\DirectedLick\B6J635_SAREmCloverRC107\20170725\1\F_B6J635_SAREmCloverRC107_20170725_plane1_Ch1_Nk26_proc.mat';
% file{2}='E:\home\ImagingData\DirectedLick\B6J635_SAREmCloverRC107\20170801\1_3\F_B6J635_SAREmCloverRC107_20170801_plane1_Ch1_Nk26_proc.mat';
% file{3}='E:\home\ImagingData\DirectedLick\B6J635_SAREmCloverRC107\20170802\1_2_3\F_B6J635_SAREmCloverRC107_20170802_plane1_Ch1_Nk26_proc.mat';

CaRoiGu_v02; % load a db, scan it.
h=guidata(CaRoiGu_v02);
%%
Nentry = length(h.DB.db)
ProcSpkDone = zeros(1,Nentry);
for ii=1:Nentry
    ProcSpkDone(ii)=h.DB.db(ii).scan_results.process.ProcSpkDone;
end
ProcSpkDone_index = find(ProcSpkDone);

% MultiIMG = [];
dats = [];
cnt=1;
for ii=1:length(ProcSpkDone_index)
    ProcSpkPath = h.DB.db(ProcSpkDone_index(ii)).scan_results.ResultsSavePath;
    ProcSpkFile = h.DB.db(ProcSpkDone_index(ii)).scan_results.ProcSpkFile;
    FullProcSpkFile = fullfile(ProcSpkPath,ProcSpkFile);
    fprintf('Loading %s...\n',FullProcSpkFile)
    dat=load(FullProcSpkFile,'cl','res','stat','ops');
    dats{ii}=dat;
    iscell=dat.cl.selected;
        
    dats{ii}.res.selected_roimask   = iscell(dat.res.iclust1);
    dats{ii}.res.roi_localind   =  dats{ii}.res.iclust1.* dats{ii}.res.selected_roimask;
    
    
end
    
%% plot them
myfigure('Mean Img');
set(gcf,'Position',[490 200 990 500]);
for ii=1:length(dats)
    subplot(2,length(dats),ii);
    imagesc(dats{ii}.res.M0);
end

for ii=1:length(dats)
    subplot(2,length(dats),length(dats)+ii);
    imagesc(dats{ii}.res.selected_roimask);
end

%%
myfigure('Overlaid');
set(gcf,'Position',[490 200 990 500]);

[optimizer,metric] = imregconfig('multimodal');
optimizer.InitialRadius = optimizer.InitialRadius/3.5;
optimizer.MaximumIterations = 300;
Fixed = dats{1}.res.M0;
moved_target = [];

% First, get the transformation 
for ii=2:length(dats)
    fprintf('%d/%d...',ii,length(dats))
    Target = dats{ii}.res.M0;
    dats{ii}.tform = imregtform(Target, Fixed, 'affine', optimizer, metric);
end
 fprintf('Done\n')
% tmp= zeros(size(dats{1}.res.roi,1), size(dats{1}.res.roi,2),3);
% tmp(:,:,1)=256*edge(dats{1}.res.roi);
%% register by mean image
dats{1}.res.roi_localind_reg=dats{1}.res.roi_localind;
dats{1}.res.M0_reg         = dats{1}.res.M0;
 dats{1}.res.probabilities_reg = dats{1}.res.probabilities;
for ii=2:length(dats)
    subplot(1,length(dats)-1,ii-1);
    Target = dats{ii}.res.M0;
    Target_moved = imwarp(Target,dats{ii}.tform,'OutputView',imref2d(size(Fixed)));
    dats{ii}.res.M0_reg = Target_moved;
    dats{ii}.res.roi_localind_reg = imwarp(dats{ii}.res.roi_localind, dats{ii}.tform,'nearest',...
        'OutputView',imref2d(size(Fixed)));
    dats{ii}.res.probabilities_reg =imwarp(reshape(dats{ii}.res.probabilities,dats{ii}.res.Ly,dats{ii}.res.Lx), dats{ii}.tform,'nearest',...
        'OutputView',imref2d(size(Fixed)));
    imshowpair(Target_moved, Fixed,'falsecolor','ColorChannels','red-cyan');
    title('Cyan: fixed, Red: moved target');
%     drawnow;
%     subplot(2,length(dats)-1,length(dats)+ii-2);
%     Target = MultiIMG(ii).ROIs;
%     Target_moved = imwarp(Target,dats{ii}.tform,'OutputView',imref2d(size(Fixed)));
%     tmp(:,:,2)=256*edge(Target_moved);
%     imshow(tmp);
    
end
%% calculate boundaries of each ROI 
for ii=1:length(dats)
    clf;
    ROI_IDs = sort(unique(dats{ii}.res.roi_localind_reg));
    ROI_IDs(ROI_IDs==0)= []; % remove 0
    ROI_simplified = zeros(size(dats{ii}.res.roi_localind_reg));
    
    se=strel('disk',3,4);
 
    B = [];
    CoM = [];
    % simplify the ROI region.
    for jj=1:length(ROI_IDs)
        
        tmp =   dats{ii}.res.roi_localind_reg==ROI_IDs(jj);
        [ind_y,ind_x]=find(tmp);
        
        CoM{jj} = round(mean([ind_y ind_x],1));
        B{jj}=bwtraceboundary(tmp,[ind_y(1), ind_x(1)],'N');
        ROI_simplified = ROI_simplified+jj*tmp;
    end

    dats{ii}.res.B = B;
    dats{ii}.res.ROI_simplified = ROI_simplified;
    dats{ii}.res.CoM = CoM;
    %     imagesc(ROI_simplified);  hold on;
    %     for jj=1:length(ROI_IDs)
    %         plot(B{jj}(:,2),B{jj}(:,1),'r'); hold on;
    %     end
    %     title(sprintf('Session %d',ii));
    %     pause(0.1);
end
  %% Then, overlay plot all the boundaries.
md_roi_h = myfigure('MultiDay ROIs');    clf; 
%     BackGroundIMG = MultiIMG(3).moved_mimg;
 BackGroundIMG = dats{ii}.res.M0_reg;
    imagesc(BackGroundIMG);colormap gray; 
    hold on;
    
    ColorOrder = get(gca,'ColorOrder');
%     LineWidthOrder = [1 1 1];
for ii=1:length(dats)
    for jj=1:length(dats{ii}.res.B)
        plot(dats{ii}.res.B{jj}(:,2),dats{ii}.res.B{jj}(:,1),...
            'Color',ColorOrder(ii,:),'LineWidth',1); hold on;
%         plot(MultiIMG(ii).CoM{jj}(2),MultiIMG(ii).CoM{jj}(1),'o','Color',ColorOrder(ii,:)); hold on;
    end  
end
set(gca,'YDir','reverse');
mean_pixel_F = mean(BackGroundIMG(:));
std_pixel_F  = std(BackGroundIMG(:));
Clim_range = mean_pixel_F+[-1 3]*std_pixel_F;
Clim_range(1)=max(0,Clim_range(1));
caxis(Clim_range);

%% Give a global ID for each ROI of each day.

% Algorithm: 
% 1. initialization.  For a given ROI plane (reference), give a global ID. 
% 2. label global ID. For the remaining ROI planes, calculate overwrap with reference 
% and give the global ID, if overlap > threshold 
% 3. delete global rois. For the remaining ROI planes, remove the ROIs that have global ID.
% This procedure leaves orphan ROIs and removes the ROIs that exists over multiple sessions.
% 
% Exclude the reference plane from ROI plane set, and repeat 1-2-3 until no
% plane remains.


Ns = length(dats);
NglobalID= 0;
% ID_table: rows are global ID (1st row is roi#1 in global, 2nd row is
% roi#2, ...). columns are sessions. 
% values in the matrix is local_ID.
% For example, ID_table(n,d) means global_ID <n> has the local ID <ID_table(n,d)> on session <d>.
ID_table = []; 
% init 
Overlap_threshold = 0.3;
for ii=1:Ns
    dats{ii}.res.roi_tmp = dats{ii}.res.roi_localind_reg;
    dats{ii}.res.roi_global_id = zeros(dats{ii}.res.Ly,dats{ii}.res.Lx);
end
myfigure('Multiday ROI');clf;
cnt = 1;
for ii=1:Ns
    
    ref_img = dats{ii}.res.roi_tmp;
    ref_roi_labels = unique(ref_img);ref_roi_labels(ref_roi_labels==0)=[];
    nROI_reg = length(ref_roi_labels);
    
    globalID = NglobalID+[1:length(ref_roi_labels)];
    
     for jj=ii:Ns
         axH(cnt)=mysubplot(Ns,Ns,ii,jj);cla;
              
         target_img = dats{jj}.res.roi_tmp;
         target_roi_labels = unique(target_img);target_roi_labels(target_roi_labels==0)=[];
         nROI_tgt = length(target_roi_labels);
         
         [Overlap,Iref,Itgt] = ROIOverlapMatrix_001(ref_img,target_img);
         
         % find roi in target that is maximally overlapped with roi in ref.
         [max_c,max_ind]=max(Overlap,[],2);
         
         % put localID in ID_table 
         idtmp = zeros(length(max_c),1);
         overlap_ok= max_c>Overlap_threshold;
         idtmp(overlap_ok)=Itgt(max_ind(overlap_ok));
         
         ID_table(globalID,jj)=idtmp;
         
         matched_id = sort(Itgt(max_ind(overlap_ok)));
         unmatched_id = setdiff(Itgt,matched_id);
         
         % remove matched id
         opt.hue_range = [0.3 0.31];
         opt.highlight_labels= matched_id;
         opt.highlight_color = 1;
         I = ROI_gem_img(dats{jj}.res.roi_localind_reg,dats{jj}.res.probabilities_reg,dats{jj}.res.M0_reg,opt);
         imagesc(I);
         axis off; 
         if ii==1, title(dats{jj}.ops.date,'Interpreter','none'), end
         
         dats{jj}.res.roi_tmp=swap_labels_IMG(dats{jj}.res.roi_tmp,matched_id,zeros(1,length(matched_id)));
            
            cnt = cnt+1;
     end
     NglobalID = length(globalID);
end

linkaxes(axH,'xy')
%%
ops1 = h.DB.ops;
savepath = fullfile(ops1.RootStorage,ops1.mouse_name,'Summary');
if ~exist(savepath,'dir'), mkdir(savepath); end

savename = sprintf('%s_multisession_reg',ops1.mouse_name);
% mysubplot(Ns,Ns,Ns,1);cla;
text(-(Ns-1)*diff(xlim),mean(ylim),savename,'interpreter','none');

print(gcf,'-dpng',fullfile(savepath,[savename,'.png']));
print(gcf,'-djpeg','-r300',fullfile(savepath,[savename,'.jpg']));
