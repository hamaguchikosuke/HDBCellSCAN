function RegImg = Load_register_OtherChannels_20181001(h)
% merge CTB and red-beads data onto tiff movie
% load('N:\home\ImagingData\DualLickMice\B6N1189_GC6s_redbPons_CTB647hM4DiRALM\20180921_d32\1_592um\ROI_B6N1189_GC6s_redbPons_CTB647hM4DiRALM_20180921_d32_plane1_Ch1_MinClust95.mat')

% load('N:\home\ImagingData\DualLickMice\B6N1189_GC6s_redbPons_CTB647hM4DiRALM\20180927_d35\1_555um\ROI_B6N1189_GC6s_redbPons_CTB647hM4DiRALM_20180927_d35_plane1_Ch1_MinClust95.mat');

ops= h.dat.ops;

[LoadPath,LoadName,Ext]=fileparts(ops.ProcFileName);
FullRGBFile = fullfile(LoadPath,[LoadName,'_TracerOverlaid.mat']); 
if exist(FullRGBFile,'file')
    answer=questdlg(sprintf('Previously registered image found (%s). Load it?',FullRGBFile));
   switch answer
       case 'Yes'
           fprintf('User selected to load previously registered data\n');
               tmp=load(FullRGBFile,'RegImg');
               RegImg=tmp.RegImg;
               return;
       otherwise
            fprintf('User selected to register images from scratch\n');
   end
end
    
%% 
MainPath =ops.RegTiffPath{1}; 
PlaneChPath = sprintf('plane%d_ch%d',ops.PlaneID,ops.ChannelID);
TiffPath = fullfile(MainPath,PlaneChPath,'x4movie');
S=dir(fullfile(TiffPath,'*.tif'));
FullTiffMoviePath = fullfile(TiffPath,S(1).name);


% TiffMovie = 'plane1_ch1\x4movie\20180921_d32_1_592um_B6N1189_GC6s_redbPons_CTB647hM4DiRALM_2P_plane1_ch1_x4_005.tif';
% try
%     RCy5Path = fullfile(MainPath,'CTB_redbeads');
%     RCy5File    = 'AVG_B6N1189_20180921_830nm_FRCy5Max26.1x[0.5]dBm_A2.5L1.5[Ave].tif';
% catch
[RCy5File,RCy5Path]=uigetfile(fullfile(MainPath,'*.*'),'Please select the file that contains red and far-red');
% endh
% try
%     GRPath = fullfile(MainPath,'CTB_redbeads');
%     GRFile    = 'AVG_B6N1189_20180921_830nm_FGRMax26.1x[0.4]dBm_A2.5L1.5[Ave].tif';
% catch
[GRFile,GRPath]=uigetfile(fullfile(RCy5Path,'*.*'),'Please select the file that contains green and red');
% end
    
meanF =h.dat.graph.mimg(:,:,2);
RCy5F= BigTiffReader(fullfile(RCy5Path,RCy5File));
GRF = BigTiffReader(fullfile(GRPath,GRFile));
%% if RCy5F and GRF is 3D image,
% [W,H,Z]=size(RCy5F);
% if Z>1
%     Nch=2;
%     RCy5F=reshape(RCy5F,W,H,2,[]);
% end
% [dsnew, Corr]  = registration_offsets(RCy5F, ops, 1);
% RCy5F  = register_movie(RCy5F, ops, dsnew);
% 
% [dsnew, Corr]  = registration_offsets(GRF, ops, 1);
% GRF  = register_movie(GRF, ops, dsnew);

%% now align to the mean F (reference and target) with Green,
% then register Red using Green-Red combination image. 
% And then, align Red-FarRed image to Red channel of Green-Red image. 

[W,H,Nch1]=size(GRF);
AllFigsH=[];
figh=myfigure('Green-Red');clf;AllFigsH=[AllFigsH,figh];
screensize=get(0,'ScreenSize');
WH=1000;WH=min([WH,screensize(3:4)]);
set(figh,'Position',[0 screensize(4)-WH WH WH]);
subplot(2,Nch1,1);
imshow(meanF);title('Fixed target (Ca)');
caxis('auto');

for ii=1:Nch1
subplot(2,Nch1,Nch1+ii); 
imshow(GRF(:,:,ii));title(sprintf('(%d)',ii));
caxis('auto')
end
colormap gray
%%
prompt = {'Align Ca data first. Select Ca image channel.'};
dlg_title = 'Input';
num_lines = 1;
defaultans = {'1'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
selected_Ch=str2double(answer{1});
Unselected_Ch = setdiff([1:Nch1],selected_Ch);


[optimizer,metric] = imregconfig('multimodal');
optimizer.InitialRadius = optimizer.InitialRadius;
optimizer.MaximumIterations = 300;
Fixed = meanF;
Target = GRF(:,:,selected_Ch);
clf;
imshowpair(Target, Fixed,'falsecolor','ColorChannels','red-cyan');
title('Before');

%
fprintf('Linear translation to match...\n')
tform_FGR = imregtform(Target, Fixed, 'translation', optimizer, metric);
fprintf('done\n');

Target_moved = imwarp(Target,tform_FGR,'OutputView',imref2d(size(Fixed)));
imshowpair(Target_moved, Fixed,'falsecolor','ColorChannels','red-cyan');
title('Cyan: fixed, Red: moved target');

GRF_moved=zeros([size(Fixed),Nch1]);

for ii=1:size(GRF,3)
    GRF_moved(:,:,ii) = imwarp(GRF(:,:,ii),tform_FGR,'OutputView',imref2d(size(Fixed)));
end
%%

[W,H,Nch2]=size(RCy5F);

figh=myfigure('Red-Cy5');clf;AllFigsH=[AllFigsH,figh];
set(figh,'Position',[0 screensize(4)-WH WH WH]);

Fixed_and_moved = cat(3,meanF,GRF_moved(:,:,Unselected_Ch));

subplot(2,Nch2,1);
imshow(Fixed_and_moved(:,:,1));title('Fixed Target 1');
caxis('auto');

subplot(2,Nch2,2);
imshow(Fixed_and_moved(:,:,2));title(sprintf('Fixed Target2: Ch%d in GR',Unselected_Ch));
caxis('auto');


for ii=1:Nch2
subplot(2,Nch2,Nch2+ii); 
imshow(RCy5F(:,:,ii));title(sprintf('Moving(%d)',ii));
caxis('auto')
end
colormap gray
%% align Red-FarRed image to Red channel of Green-Red image. 

prompt = {'Which channel to use as a target from "Fixed Target"?','Which channel to use as Moving channel?'};
dlg_title = 'Input';
num_lines = 1;
defaultans = {'1','2'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

selected_FixedCh=str2double(answer{1});
Unselected_FixedCh = setdiff([1:Nch1],selected_FixedCh);

selected_MovingCh=str2double(answer{2});
Unselected_MovingCh = setdiff([1:Nch1],selected_MovingCh);

Fixed = Fixed_and_moved(:,:,selected_FixedCh);
Moving = RCy5F(:,:,selected_MovingCh);
clf;
imshowpair(Moving, Fixed,'falsecolor','ColorChannels','red-cyan');
title('Before');
drawnow;

fprintf('Linear translation to match...\n')
tform_RCy5 = imregtform(Moving, Fixed, 'translation', optimizer, metric);
fprintf('done\n');



RCy5_moved=zeros([size(Fixed),Nch2]);

for ii=1:size(GRF,3)
    RCy5_moved(:,:,ii) = imwarp(RCy5F(:,:,ii),tform_RCy5,'OutputView',imref2d(size(Fixed)));
end

Target_moved=RCy5_moved(:,:,selected_MovingCh);
imshowpair(Target_moved, Fixed,'falsecolor','ColorChannels','red-cyan');
title('Done');
%%
% figh=myfigure('GC6s_redbeads_CTB_movie');clf;

% set(gcf,'Position',[50 500 1050 200]);

figh=myfigure('RCy5 and GRF');clf;
set(figh,'Position',[0 screensize(4)-WH/2 2*WH WH/2]);AllFigsH=[AllFigsH,figh];
RGCy5Img_moved = cat(3,GRF_moved,RCy5_moved);

for ii=1:size(RGCy5Img_moved,3)
subplot(1,size(RGCy5Img_moved,3),ii); 
imshow(RGCy5Img_moved(:,:,ii));title(sprintf('(%d)',ii));
caxis('auto');
end
colormap gray

prompt = {'Select red channel data to overlay with movie', 'blue channel'};
dlg_title = 'Input';
num_lines = 1;
defaultans = {'1','2'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

red_ch  =str2double(answer{1});
blue_ch =str2double(answer{2});
%%
R=RGCy5Img_moved(:,:,red_ch);
B=RGCy5Img_moved(:,:,blue_ch);
G=meanF;

r_max = max(R(:));
r_bias =min(R(:));

g_max = max(G(:));
g_bias = min(G(:));

b_max = max(B(:));
b_bias = min(B(:));

% maxRGB=max([r_max,g_max,b_max]);
r_coef = 1/r_max;
g_coef= 1/g_max;
b_coef = 1/b_max;


Composite=cat(3,...
    r_coef*(R-r_bias),...
    g_coef*(G-g_bias),...
    b_coef*(B-b_bias));
clf;
imshow(Composite);
title(sprintf('Overlaid %s',FullTiffMoviePath),'Interpreter','none','FontSize',7);
[SavePath,SaveName,Ext]=fileparts(ops.ProcFileName);

 FullTiffFile = fullfile(SavePath,[SaveName,'_TracerOverlaid.tif']);
 fprintf('Printing %s...\n',FullTiffFile);
 print(gcf,'-dtiff',FullTiffFile);
 

 
 RegImg.R = R;
 RegImg.G = G;
 RegImg.B = B;
 RegImg.Ca= 2;
 RegImg.Composite=Composite;
 
 FullMATFile = fullfile(SavePath,[SaveName,'_TracerOverlaid.mat']);
 save(FullMATFile,'RegImg');
 
 
 delete(AllFigsH);
%  TiffWriter(Composite,FullTiffFile,16);
%% BackG = RedBeadsCTB;
% clf;
% movieF=zeros(W,H,3,T,'uint16');
% FullVideoFile = fullfile(MainPath,'TracerOverlaied.avi');
% 
% myVideo = VideoWriter(FullVideoFile);
% myVideo.FrameRate = 15;  % Default 30
% myVideo.Quality = 100;    % Default 75
% open(myVideo);
% for ii=1:T
%     Composite=cat(3,r_coef*(BackG(:,:,red_ch)+r_bias),...
%         g_coef*(gF(:,:,ii)+g_bias),...
%         b_coef*(BackG(:,:,blue_ch))+b_bias);
%     imshow(Composite);title(sprintf('frame=%d',ii));
%     pause(0.01);
%     frame = getframe(gca);
%     writeVideo(myVideo,frame);
% end
% 
% close(myVideo);
% 

