%% align ROIs from different days 
% 1. load pipeline-processed proc data
%    load the mean data 
%  <To Do>: change this to proc-behav data.
file{1}='E:\home\ImagingData\DirectedLick\B6J635_SAREmCloverRC107\20170725\1\F_B6J635_SAREmCloverRC107_20170725_plane1_Ch1_Nk26_proc.mat';
file{2}='E:\home\ImagingData\DirectedLick\B6J635_SAREmCloverRC107\20170801\1_3\F_B6J635_SAREmCloverRC107_20170801_plane1_Ch1_Nk26_proc.mat';
file{3}='E:\home\ImagingData\DirectedLick\B6J635_SAREmCloverRC107\20170802\1_2_3\F_B6J635_SAREmCloverRC107_20170802_plane1_Ch1_Nk26_proc.mat';


MultiIMG = [];
dats = [];
for ii=1:length(file)
    fprintf('Loading %s...\n',file{ii})
    load(file{ii},'dat');
    dats{ii}=dat;
    iscell=dat.cl.iscell;
    iscell_numbered = iscell;
    iscell_numbered(find(iscell_numbered))=1:length(find(iscell_numbered));
    
    MultiIMG(ii).ROIs          =reshape(iscell(dat.res.iclust),          dat.cl.Ly, dat.cl.Lx);
    MultiIMG(ii).ROIs_numbered =reshape(iscell_numbered(dat.res.iclust), dat.cl.Ly, dat.cl.Lx);
    MultiIMG(ii).ming =dat.mimg;
end
    
%% plot them
myfigure('Mean Img');
set(gcf,'Position',[490 200 990 500]);
for ii=1:length(dats)
    subplot(2,length(dats),ii);
    imagesc(dats{ii}.mimg(:,:,2));
end

for ii=1:length(dats)
    subplot(2,length(dats),length(dats)+ii);
    imagesc(MultiIMG(ii).ROIs);
end

%%
myfigure('Overlaid');
set(gcf,'Position',[490 200 990 500]);

[optimizer,metric] = imregconfig('multimodal');
optimizer.InitialRadius = optimizer.InitialRadius/3.5;
optimizer.MaximumIterations = 300;
Fixed = dats{1}.mimg(:,:,2);
moved_target = [];

% First, get the transformation 
for ii=2:length(dats)
    Target = dats{ii}.mimg(:,:,2);
    dats{ii}.tform = imregtform(Target, Fixed, 'affine', optimizer, metric);
end

tmp= zeros(size(MultiIMG(1).ROIs,1), size(MultiIMG(1).ROIs,2),3);
tmp(:,:,1)=256*edge(MultiIMG(1).ROIs);
%% register by mean image
 MultiIMG(1).ROIs_numbered_moved=MultiIMG(1).ROIs_numbered;
 MultiIMG(1).moved_mimg = dats{1}.mimg(:,:,2);;
 
for ii=2:length(dats)
    subplot(1,length(dats)-1,ii-1);
    Target = dats{ii}.mimg(:,:,2);
    Target_moved = imwarp(Target,dats{ii}.tform,'OutputView',imref2d(size(Fixed)));
    MultiIMG(ii).moved_mimg = Target_moved;
    MultiIMG(ii).ROIs_numbered_moved = imwarp(MultiIMG(ii).ROIs_numbered, dats{ii}.tform,'nearest',...
        'OutputView',imref2d(size(Fixed)));
    imshowpair(Target_moved, Fixed,'falsecolor','ColorChannels','red-cyan');
    title('Cyan: fixed, Red: moved target');
    
%     subplot(2,length(dats)-1,length(dats)+ii-2);
%     Target = MultiIMG(ii).ROIs;
%     Target_moved = imwarp(Target,dats{ii}.tform,'OutputView',imref2d(size(Fixed)));
%     tmp(:,:,2)=256*edge(Target_moved);
%     imshow(tmp);
    
end
%% calculate boundaries of each ROI 
for ii=1:length(MultiIMG);
    clf;
    ROI_IDs = sort(unique(MultiIMG(ii).ROIs_numbered_moved));
    ROI_IDs(ROI_IDs==0)= []; % remove 0
    ROI_simplified = zeros(size(MultiIMG(ii).ROIs_numbered_moved));
    
    se=strel('disk',3,4);
 
    B = [];
    CoM = [];
    % simplify the ROI region. 
    for jj=1:length(ROI_IDs)
        %         jj
        ind_x=[];
        ind_y = [];
        init_disk_size = 5;
        % in case the ROI is too small and imerode eliminate it, reduce the
        % size of erosion. 
        while isempty(ind_y)
            init_disk_size = init_disk_size-1;
            se=strel('disk',init_disk_size,4);
            tmp = MultiIMG(ii).ROIs_numbered_moved==ROI_IDs(jj);
            tmp = imerode(tmp,se,'same');
            tmp = imdilate(tmp,se,'same');
            [ind_y,ind_x]=find(tmp);
            if init_disk_size==0
                tmp =   MultiIMG(ii).ROIs_numbered_moved==ROI_IDs(jj);
                 [ind_y,ind_x]=find(tmp);
                return;
            end
        end
        CoM{jj} = round(mean([ind_y ind_x],1));
        B{jj}=bwtraceboundary(tmp,[ind_y(1), ind_x(1)],'N');
        %     plot(B{jj}(:,2),B{jj}(:,1),'r'); hold on;
        ROI_simplified = ROI_simplified+jj*tmp;
    end
    % imagesc(MultiIMG(ii).ROIs_numbered);
        imagesc(ROI_simplified);  hold on;
        for jj=1:length(ROI_IDs)
            plot(B{jj}(:,2),B{jj}(:,1),'r'); hold on;
        end
        MultiIMG(ii).B = B;
        MultiIMG(ii).ROI_simplified = ROI_simplified;
        MultiIMG(ii).CoM = CoM;
        pause(1);
end
  %% Then, overlay plot all the boundaries.
md_roi_h = myfigure('MultiDay ROIs');    clf; 
    BackGroundIMG = MultiIMG(3).moved_mimg;
    imagesc(BackGroundIMG);colormap gray; 
    hold on;
    
    ColorOrder = get(gca,'ColorOrder');
for ii=1:length(MultiIMG);
    for jj=1:length(MultiIMG(ii).B)
        plot(MultiIMG(ii).B{jj}(:,2),MultiIMG(ii).B{jj}(:,1),'Color',ColorOrder(ii,:),'LineWidth',1.5); hold on;
%         plot(MultiIMG(ii).CoM{jj}(2),MultiIMG(ii).CoM{jj}(1),'o','Color',ColorOrder(ii,:)); hold on;
    end  
end
set(gca,'YDir','reverse');
mean_pixel_F = mean(BackGroundIMG(:));
std_pixel_F  = std(BackGroundIMG(:));
Clim_range = mean_pixel_F+[-1 3]*std_pixel_F;
Clim_range(1)=max(0,Clim_range(1));
caxis(Clim_range);
%% for a given pixel point, find the ROI that contains that pixel.
figure(md_roi_h);
[ind_x,ind_y]=ginput(1);
ind_x = round(ind_x);
ind_y = round(ind_y);
fprintf('Targeted [%d,%d]\n',ind_x,ind_y);

myfigure('F Plot over days');clf;
for ii=1:length(MultiIMG)
    ROI_id = MultiIMG(ii).ROIs_numbered_moved(ind_y,ind_x);
    fprintf('%d: ROI %d \n',ii,ROI_id);
    plot()
% MultiIMG(ii).ROIs_numbered_moved(ind_y,ind_x)
end
%% now, I can give a global ID for each ROI of each day.
% for example, 
% dats{ii}.ROI(jj)

pos=get(gca,'CurrentPoint');

x = round(pos(1,1));
y  = round(pos(1,2));

x = min(max(1, round(x)), h.dat.cl.Lx);
y = min(max(1, round(y)), h.dat.cl.Ly);

h.dat.F.ichosen = h.dat.res.iclust(y, x);


%% additional fun stuff.

hFig=myfigure('test');
set(hFig,'Toolbar','none','Menubar','none');
hIm = imshow(tmp);
hSP = imscrollpanel(hFig,hIm);
set(hSP,'Units','normalized',...
	  'Position',[0 .1 1 .9])

hMagBox = immagbox(hFig,hIm);
pos = get(hMagBox,'Position');
set(hMagBox,'Position',[0 0 pos(3) pos(4)])
imoverview(hIm);


% 3. Get the scroll panel API to programmatically control the view.
api = iptgetapi(hSP);

% 4. Get the current magnification and position.
mag = api.getMagnification();
r = api.getVisibleImageRect();

% 5. View the top left corner of the image.
api.setVisibleLocation(0.5,0.5)

% 6. Change the magnification to the value that just fits.
api.setMagnification(api.findFitMag())

% 7. Zoom in to 1600% on the dark spot.
api.setMagnificationAndCenter(16,306,800)


