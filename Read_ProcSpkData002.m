% function out = Read_ProcSpkData002()
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

% AnimalIDs=[3250,3479,3482];opt='all'

%% RdG High 
% AnimalIDs=[3241,3249,3250,3479,3482]; opt='all';suffix='';GroupName='High(jG8m)';
%% RdG Low 
% AnimalIDs=[3242,3251,3252,3572]; opt='all';suffix='';GroupName='Low(jG8m)'; % selected the plane that was imaged through the whole 16 weeks.

%% EnvB RV->training 
AnimalIDs=[3344,3346,3420,3506];opt='spont';suffix='_spont';GroupName='EnvB(RV-3w-Training)';

%% EnvB Training&RV
% AnimalIDs=[3414,3512,3519,3520,3615,3616];opt='spont';suffix='_spont';GroupName='EnvB(Training_and_RV)';
% AnimalIDs=[3250];
% AnimalIDs=[3414];opt='spont';


xgtb=[];
for aa=1:length(AnimalIDs)      
    AnimalID=AnimalIDs(aa);
    Get_SamePlaneFcn=eval(sprintf('@get_SamePlane_%d%s',AnimalID,suffix));
    SamePlane=get_SamePlaneWrapperFcn(Get_SamePlaneFcn,opt);

    for ii=1:length(SamePlane)
        FileName=SamePlane(ii).FullProcSpkFileS;
        
     
        if isempty(FileName), continue; end
        if length(FileName)==1, continue; end
        [x,w,xlabelname]=localanalysis(FileName);
        xgt=array2table(x,'VariableNames',{'activity'});
        xgt.week=w;
        xgt.animal=AnimalID*ones(height(xgt),1);
        xgt.plane=ii*ones(height(xgt),1);
        xgtb=cat(1,xgtb,xgt);
    end
end

weekvalueset = string([1:16]);
weekvalueset=[weekvalueset,"noise"];
xgtb.week=categorical(xgtb.week,weekvalueset,Ordinal=true);

s=strsplit(FileName{1},filesep);
SaveDir=fullfile(s{1:end-4},'Summary');
mkdir(SaveDir);
FileName=sprintf('%d_',AnimalIDs);
FileName=FileName(1:end-1);
fullFileName=fullfile(SaveDir,sprintf('%s_%s',GroupName,FileName));

writetable(xgtb,[fullFileName,'.csv']);

%%

myfigure('Activity');
setPaperSize([8,12]);clf;
CountsPerPlane=grpstats(xgtb,{'week','animal','plane'},{'numel'});

mxgtb=grpstats(xgtb,{'week','animal'},{'mean','sem'},'DataVars',{'activity'});
mxgtb=sortrows(mxgtb,{'animal','week'},{'ascend','ascend'});


uniqMouseIDs=unique(mxgtb.animal);
Nmouse=length(uniqMouseIDs);
h=linspace(0,1,Nmouse+1)';
h=h(1:end-1);
rng(20241206);
h=h(randperm(length(h)));
hsv=[h,ones(Nmouse,2)];
rgb=hsv2rgb(hsv);
colo=[];
colo.animal=uniqMouseIDs;
colo.rgb=rgb;
colo=struct2table(colo);
mxgtb=innerjoin(mxgtb,colo,keys={'animal'})

subH=[];
for animalID=AnimalIDs
    subH(1)=subplot(2,1,1);
    data=mxgtb(mxgtb.animal==animalID,:);
    h=plot(data.week,data.mean_activity,'o'); hold on;
    set(h,'Color',data.rgb(1,:));
    %     errorbar(data.week,data.mean_x,data.Fun3_x,data.Fun4_x);
    h=errorbar(data.week(1:end-1),data.mean_activity(1:end-1),3*data.sem_activity(1:end-1));
    set(h,'Color',data.rgb(1,:));
    h=errorbar(data.week(end),data.mean_activity(end),3*data.sem_activity(end));
    set(h,'Color',data.rgb(1,:));


    ylabel('Mean Activity');
    title([GroupName,':',sprintf('%d_',AnimalIDs)],'Interpreter','none')
    subH(2)=subplot(2,1,2);
    h=plot(data.week(1:end-1),data.GroupCount(1:end-1),'o-'); hold on; % no need t
    h.UserData=data;
%     report_userdata=@(src,~)disp(src.UserData);
    h.ButtonDownFcn=@report_userdata;
    
    set(h,'Color',data.rgb(1,:));
    set(gca,'YLim',[0 250]);
    ylabel('# ROI detected')
    xlabel('week');
end
set(subH,'XTickLabelRotation',0);

fprintf('Printing %s....\n',[fullFileName,'_population_activity'])
print(gcf,'-dpdf',[fullFileName,'_population_activity.pdf'])
print(gcf,'-dpng',[fullFileName,'_population_activity.png'])


function report_userdata(src,~)
    disp(src.UserData);
end

%%
function [x,w,xlabelname]=localanalysis(FileName)
% --input--
% FileName: a cell of FullFileNames designating the proc.mat file
% --output--
% xg: values and group labels.
% xg(:,1) contains values 
% xg(:,2) contains labels (integers)
% xlabelname is the name of the groups 
% 

N=length(FileName);
xlabelname={};
for ii=1:length(FileName)
    s=strsplit(FileName{ii},filesep);
    xlabelname{ii}=strrep(s{end-2},'_','\_');
end
xlabelname{end+1}='noise';
animLabel=s{end-3};
planeLabel=s{end-1};

weekcnt=regexp(xlabelname,'.*week(?<week>\d+)','names');
if isempty(weekcnt{1})
    weekcnt=regexp(xlabelname,'.*_(?<week>\d+)w','names');
end
weekcnt=[weekcnt{:}];
weekcnt={weekcnt.week,'noise'};


%%
% MeanActivity1w={};
MeanActivity={};
MeanActivityNoise={};

zx={};
ts={};
mimg={};

Labels=[];
for ii=1:N
    s=strsplit(FileName{ii},filesep);
    dateLabel=s{end-2};
    Labels{ii}=dateLabel;
end
SaveDir=fullfile(s{1:end-3},'Summary');

Fspk='F';
for ii=1:N
    Data  = GetMeanActivity(FileName{ii},Fspk);
    MeanActivity{ii}=Data.summary.MeanActivity;
    MeanActivityNoise{ii}=Data.summary.MeanActivityNoise;
%     Fs{ii}=Data.summary.F;
%     spks{ii}=Data.summary.spk;
    zx{ii}=Data.summary.zx(Data.cl.selected>0,:);
    ts{ii}=Data.summary.t;
    mimg{ii}=Data.graph.mimg(:,:,2);
end
% %%
% histH=[];
% figure(1);clf;
% % colors={'b','g','r'};
% 
% for ii=1:N
% %     histH(ii)=histogram(MeanActivity{ii},linspace(0,0.5,30),'Normalization','probability');hold on; set(histH(ii),'FaceColor',colors{ii});
%     histH(ii)=histogram(MeanActivity{ii},linspace(0,0.5,30),'Normalization','probability');hold on; set(histH(ii));
% end
% legend(histH,{Labels{1:N}});
% xlabel('Mean Activity [time average of variance-normalized deconvoluted F - its median]');
% ylabel('Prob.')
%% Example of Traces that have high Mean Activity 
myfigure('stackedPlot');
setPaperSize([6,24]);clf;
yshift=0;
shiftsize=5;

GoodThreshold = 0.25;
GoodThresholdMax = inf;
hsv=linspace(0,1,N+1);hsv=hsv(1:end-1);
colors=hsv2rgb([hsv(:),ones(N,2)]);

for ii=1:N
    GoodOnes=MeanActivity{ii}>GoodThreshold & MeanActivity{ii}<GoodThresholdMax;
    if any(GoodOnes)
    xplots = zx{ii}(GoodOnes,:);
    xplots = xplots+yshift-shiftsize*[1:nnz(GoodOnes)]';
    plot(ts{ii},xplots,'Color',colors(ii,:)); hold on;
    yshift = yshift-shiftsize*nnz(GoodOnes);
    end
end
ylabel(sprintf('zscore(%s)',Fspk));xlabel('sec')
title(sprintf('stacked_activity_above%g',GoodThreshold),'Interpreter','none');
axis off;

do_print=1;
if (do_print)
    fullsavename=fullfile(SaveDir,sprintf('%s_%s_stacked_activity_above%g',animLabel,planeLabel,GoodThreshold));
    fprintf('printing %s...\n',fullsavename)
    print(gcf,'-dpdf',[fullsavename,'.pdf']);
    print(gcf,'-dpng',[fullsavename,'.png']);
end

%%
histH=[];
x=[];
g=[];
w=[];
myfigure('SwarmPlot_MeanActivity');
setPaperSize([12,8]);clf;
nsamples=[];
for ii=1:length(MeanActivity)
% %             histH(ii)=histogram(MeanActivity{ii},linspace(0,12,30),'Normalization','probability');hold on; %set(histH(1),'FaceColor','b');
% %               stairs(hist_edges{ii}(1:end-1),hist_Ngood{ii});hold on;
       x=[x;MeanActivity{ii}];
       g=[g;ii*ones(length(MeanActivity{ii}),1)];
       nsamples(ii)=length(MeanActivity{ii});
end

%         legend(histH,{'1w','4w'});
%         x=cat(1,MeanActivity,MeanActivityNoise);x=reshape(x,numel(x),[]);

%         violin(x','bw',0.02);
n=cat(1,MeanActivityNoise{:});
nsamples(end+1)=length(n);
gn=(length(MeanActivity)+1)*ones(length(n),1);

% xg=cat(2,[x;n],[g;gn]);
x=[x;n];
w=weekcnt([g;gn])';

violinplot(x,w,'MarkerSize',5);
for ii=1:length(MeanActivity)+1
    text(ii-0.3,max(ylim)*0.9,sprintf('(%d)',nsamples(ii)));
end
line(xlim,[0,0])
set(gca,'XTickLabel',xlabelname);
ylabel('Activity(mean of median subtracted F)');
if (do_print)
    fullsavename=fullfile(SaveDir,sprintf('%s_%s_activity_swplt',animLabel,planeLabel));
    fprintf('Printing %s\n',fullsavename);
    print(gcf,'-dpdf',[fullsavename,'.pdf']);
    print(gcf,'-dpng',[fullsavename,'.png']);
end
%% 

myfigure('TiledPlot');
setPaperSize([12,12]);clf;
imgh=[];
Ncol=3;
Nrow=ceil(N/Ncol);
for ii=1:N
    [J,I]=ind2sub([Nrow,Ncol],ii);
    mysubplot(Nrow,Ncol,I,J,[0.98,0.9]);
    imgh(ii)=imshow(mimg{ii});
    title(Labels{ii},'interpreter','none');
    c1c99=prctile(mimg{ii}(:),[1,99.5]);
    colormap gray
    set(gca,'CLim',c1c99);
%     axis off;
end
pimgh=get(imgh,'Parent');
linkaxes([pimgh{:}],'xy')
if (do_print)
    fullsavename=fullfile(SaveDir,sprintf('%s_%s_tiled.pdf',animLabel,planeLabel));
    fprintf('Printing %s\n',fullsavename);
    print(gcf,'-dpdf',[fullsavename])
end


end

%% 
function Data= GetMeanActivity(FileName,FSpk)
   
% FSpk='F'
Data=load(FileName);
    switch lower(FSpk)
        case lower('F')
            x=Data.F.Fcell{1}-Data.F.FcellNeu{1};
        case lower('Spk')
            x=Data.F.Spk{1};
        otherwise
            error('Unknown option %s',Fspk);
    end
    x=detrend(x',1)';


FPS=30; % Frame rate
%%
sig=ceil(0.1*FPS);
GW=gausswin(6*sig);GW=GW/sum(GW);

% F=Data.F.Ftrue{1}(goodROI,:); 
% neuropilF=Data.F.FcellNeu{1}(goodROI,:);
% rawF = Data.F.Fcell{1}(goodROI,:);
% spk  = Data.F.Spk{1};
% spk(:,end-2*FPS:end)=NaN; % neglect the last 2sec which contains artifact of deconvolution.
t=[1:size(x,2)]/FPS;

do_conv = 1;
if (do_conv)
%     F=conv2(F,GW','same');
%     neuropilF=conv2(neuropilF,GW','same');
%     rawF = conv2(rawF,GW','same');
    x=conv2(x,GW','same');
end
%%
% figure(100);clf;
% axH=[];
% axH(1)=subplot(3,1,1);cla;
% ind=randi(size(F,1),1);
% 
% plot(t,[F(ind,:);neuropilF(ind,:);rawF(ind,:)]);hold on;
% legend({'true    F = rawF-neuropilF','neuropilF','rawF'});
% title(sprintf('Neuron #%d',ind))
% axH(2)=subplot(3,1,2);cla;
% plot(t,spk(ind,:),'k-');
% 
% medianSpk=quantile(spk(ind,:),0.5); % 50% percentile as a floor of activity 
% line(xlim,medianSpk*[1 1],'LineStyle','--','Color',[1 .5 .5]);
% legend('deconvolved F \approx spike','median value')
% linkaxes(axH,'x');

%% check    
% 
% axH(3)=subplot(3,1,3);cla;
    
 zx=zscore(x,0,2);
 medianzx=median(zx,2);
 dzx=zx-medianzx;
 Tmax=t(end);

% good one is manually selected ones and at least it has positive skewness
 goodROI=(Data.cl.selected>0) & (mean(dzx,2)>0);

 edges=linspace(-0.2,0.4,21);de=edges(2)-edges(1);
 MeanActivity=sum(dzx(goodROI,:),2)/Tmax/FPS;
 [Ngood]=histcounts(MeanActivity,edges);
 MeanActivityNoise=sum(dzx(~goodROI,:),2)/Tmax/FPS;
 [Nnoise]=histcounts(MeanActivityNoise,edges);
 myfigure('SNR');clf;
 stairs([edges],[0,Ngood]);hold on;
 stairs([edges],[0,Nnoise]);

 Data.summary.MeanActivity = MeanActivity;
 Data.summary.MeanActivityNoise = MeanActivityNoise;
 Data.summary.hist_Ngood=Ngood;
 Data.summary.hist_edges=edges;
 Data.summary.x=x;
 Data.summary.zx=zx;
 Data.summary.t= t;

% 
% zspk= spk./nanstd(spk,0,2); % normalized by STD so that it always have variance 1
% zmedianspk = nanmedian(zspk,2);
% DeltaSpk=zspk-zmedianspk;
% MeanActivity = FPS*nanmean(DeltaSpk,2);
% MeanActivityNoise=sum(dzx(~goodROI,:),2)/Tmax/FPS;
% histogram(MeanActivity,20);
% xlabel('Mean Activity [time average of variance-normalized deconvoluted F - its median]');
% ylabel('counts')
% 
% Data.summary.MeanActivity = MeanActivity;
% Data.summary.F=F;
% Data.summary.spk=spk;
% Data.summary.zspk=zspk;
% Data.summary.t= t;

end

