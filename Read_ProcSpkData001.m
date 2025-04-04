% function out = Read_ProcSpkData001()
%% example code to read and plot data from .procspk.mat data 
% procspk.mat contains fluorescence signal averaged within ROIs,
% and deconvoluted fluorescence signal.
%
% Data.
%       F. # Fluorecence related data
%        Ftrue # Fluorescence signal = raw ROI average - neuropil signal 
%        FcellNeu # Neuropil signal computed from the donuts like region surrounding the ROI
%        Fcell  # raw fluorescence averaged within the ROI
%        Spk    # deconvoluted fluorescence signal used as an proxy of spikes
%       cl. # cluster realted data
%         selected % ROI selected as good. Spk are computed only for selected (good) ROIs.
%
% by Kosuke Hamaguchi 20230713

% Rewrite the FileName and path accordingly. 
clear;
out=1;
FileName=[];
% B6N3246, depth 300

% FilePath = 'DirectoryPathToWhereYourDataExist';
FilePath = 'C:\Users\hammer\Desktop\209_RVG_rC26dG1-3n7-916-jGCaMP8f';
FileName{1}=fullfile(FilePath,'B6N3246_jG8f\20230616_1week\1_A2.5_L1.5_d300\Fsig_B6N3246_jG8f_20230616_1week_plane1_ch1_MinClust79_procSpk.mat');
FileName{2}=fullfile(FilePath,'B6N3246_jG8f\20230707_4week\1_A2.5_L1.5_d300\Fsig_B6N3246_jG8f_20230707_4week_plane1_ch1_MinClust79_procSpk.mat');

% B6N3246_jG8f, depth 500
FileName{1}=fullfile(FilePath,'B6N3246_jG8f\20230616_1week\1_A2.5_L1.5_d550\Fsig_B6N3246_jG8f_20230616_1week_plane1_ch1_MinClust79_procSpk.mat');
FileName{2}=fullfile(FilePath,'B6N3246_jG8f\20230707_4week\1_A2.5_L1.5_d550\Fsig_B6N3246_jG8f_20230707_4week_plane1_ch1_MinClust79_procSpk.mat');


% B6N3249_jG8m, depth 500
FilePath = 'F:\kh018\home\ImagingData\InoueRabies\#210_RVG_rC26dG1-3n7-916-jGCaMP8m\';
FileName{1}=fullfile(FilePath,'B6N3249_jG8m\20230626_1w\1_A2.5_L1.4_d500\Fsig_B6N3249_jG8m_20230626_1w_plane1_ch1_MinClust79_procSpk.mat');
FileName{2}=fullfile(FilePath,'B6N3249_jG8m\20230717_4w\1_A2.5_L1.4_d500\Fsig_B6N3249_jG8m_20230717_4w_plane1_ch1_MinClust79_procSpk.mat');


% B6N3250_jG8m, depth 500
FilePath = 'F:\kh018\home\ImagingData\InoueRabies\#210_RVG_rC26dG1-3n7-916-jGCaMP8m\';
FileName{1}=fullfile(FilePath,'B6N3250_jG8m\20230626_1w\1_A2.6_L1.5_d500\Fsig_B6N3250_jG8m_20230626_1w_plane1_ch1_MinClust79_procSpk.mat');
FileName{2}=fullfile(FilePath,'B6N3250_jG8m\20230717_4w\1_A2.6_L1.5_d500\Fsig_B6N3250_jG8m_20230717_4w_plane1_ch1_MinClust79_procSpk.mat');

% MeanActivity1w={};
MeanActivity={};
Fs={};
spks={};
ts={};
mimg={};
Labels={'1w','4w'};

ii=1;
for ii=1:2
    Data  = GetMeanActivity(FileName{ii});
    MeanActivity{ii}=Data.summary.MeanActivity;
    Fs{ii}=Data.summary.F;
    spks{ii}=Data.summary.spk;
    ts{ii}=Data.summary.t;
    mimg{ii}=Data.graph.mimg(:,:,2);
end
%%
histH=[];
figure(1);clf;
histH(1)=histogram(MeanActivity{1},linspace(0,12,30),'Normalization','probability');hold on; set(histH(1),'FaceColor','b');
histH(2)=histogram(MeanActivity{2},linspace(0,12,30),'Normalization','probability');hold on; set(histH(2),'FaceColor','g');
legend(histH,{'1w','4w'});
xlabel('Mean Activity [time average of variance-normalized deconvoluted F - its median]');
ylabel('Prob.')
%% Example of Traces that have high Mean Activity 
figure(2);clf;
yshift=0;
shiftsize=2;

GoodThreshold = 6;
GoodThresholdMax = inf;
for ii=1:2
    GoodOnes=MeanActivity{ii}>GoodThreshold & MeanActivity{ii}<GoodThresholdMax;
    if any(GoodOnes)
    spkplots = spks{ii}(GoodOnes,:);
    spkplots = spkplots+yshift-shiftsize*[1:nnz(GoodOnes)]';
    plot(ts{ii},spkplots,'Color',get(histH(ii),'FaceColor')); hold on;
    yshift = yshift-shiftsize*nnz(GoodOnes);
    end
end
ylabel('zscore(deconv F)');xlabel('sec')

figure(3);clf;
yshift=0;
shiftsize=6;
for ii=1:2
    GoodOnes=MeanActivity{ii}>GoodThreshold & MeanActivity{ii}<GoodThresholdMax;
    if any(GoodOnes)
    Fplots = zscore(Fs{ii}(GoodOnes,:),0,2);
    Fplots = Fplots+yshift-shiftsize*[1:nnz(GoodOnes)]';
    plot(ts{ii},Fplots,'Color',get(histH(ii),'FaceColor')); hold on;
    yshift = yshift-shiftsize*nnz(GoodOnes);
    end
end
ylabel('zscore(F)');xlabel('sec')


%% 
figure(4);clf;

for ii=1:2
    subplot(1,2,ii);
    imagesc(mimg{ii});title(Labels{ii});
end




%% 
function Data= GetMeanActivity(FileName)
Data=load(FileName);
FPS=30; % Frame rate
%%
goodROI=Data.cl.selected>0;
sig=ceil(0.1*FPS);
GW=gausswin(6*sig);GW=GW/sum(GW);

F=Data.F.Ftrue{1}(goodROI,:); 
neuropilF=Data.F.FcellNeu{1}(goodROI,:);
rawF = Data.F.Fcell{1}(goodROI,:);
spk  = Data.F.Spk{1};
spk(:,end-2*FPS:end)=NaN; % neglect the last 2sec which contains artifact of deconvolution.
t=[1:size(F,2)]/FPS;

do_conv = 1;
if (do_conv)
    F=conv2(F,GW','same');
    neuropilF=conv2(neuropilF,GW','same');
    rawF = conv2(rawF,GW','same');
end
%%
figure(100);clf;
axH=[];
axH(1)=subplot(3,1,1);cla;
ind=randi(size(F,1),1);

plot(t,[F(ind,:);neuropilF(ind,:);rawF(ind,:)]);hold on;
legend({'true    F = rawF-neuropilF','neuropilF','rawF'});
title(sprintf('Neuron #%d',ind))
axH(2)=subplot(3,1,2);cla;
plot(t,spk(ind,:),'k-');

medianSpk=quantile(spk(ind,:),0.5); % 50% percentile as a floor of activity 
line(xlim,medianSpk*[1 1],'LineStyle','--','Color',[1 .5 .5]);
legend('deconvolved F \approx spike','median value')
linkaxes(axH,'x');

%% check    

axH(3)=subplot(3,1,3);cla;
zspk= spk./nanstd(spk,0,2); % normalized by STD so that it always have variance 1
zmedianspk = nanmedian(zspk,2);
DeltaSpk=zspk-zmedianspk;
MeanActivity = FPS*nanmean(DeltaSpk,2);
histogram(MeanActivity,20);
xlabel('Mean Activity [time average of variance-normalized deconvoluted F - its median]');
ylabel('counts')

Data.summary.MeanActivity = MeanActivity;
Data.summary.F=F;
Data.summary.spk=spk;
Data.summary.zspk=zspk;
Data.summary.t= t;

end

