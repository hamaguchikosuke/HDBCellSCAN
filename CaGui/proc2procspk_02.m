function proc2procspk(FullProcFile)
%% load proc and add spike convolution data
%  proc2procspk(FullProcFile);
%  FullProcFile: full path and filename of *_proc.mat
%  Output is _procSpk.mat
% 
%  by KH 20171121
%  Update: 20181101: subtract_neuropil with option, Coef='1/Regress'
%  because in some neuropil dense image, I have to maximize the
%  subtraction.


if isempty(FullProcFile)
   [ProcFile,ProcPath]= uigetfile( ...
       {'*.mat';'*.*'}, ...
        'Pick a proc file');
    [~,ProcFile,ext]=fileparts(ProcFile);
else
   [ProcPath,ProcFile,ext]=fileparts(FullProcFile );
end

FullProcFile = fullfile(ProcPath,[ProcFile,ext]);
FullProcSpkFile =fullfile(ProcPath,[ProcFile,'Spk',ext]);

if exist(FullProcFile,'file')
    fprintf('Loading Proc  File (%s)\n',FullProcFile);
    dat=load(FullProcFile);
else
    error('No proc file found');
end

if isfield(dat,'dat')
    dat=dat.dat;
end

%% do Fluorescent trace deconvolution 
 
dat.F.Spk={};

NumSubDirs=length(dat.F.Fcell);

T=zeros(1,NumSubDirs);
dF = [];
for jj=1:NumSubDirs
    
    cellindex=find(dat.cl.selected);
    if max(cellindex)>size(dat.F.Fcell{jj},1)
        error('it seems that proc file is not finalized (cellindex size mismatched with dat.F.Fcell size)');
    end
    Ff=[dat.F.Fcell{jj}(cellindex,:)];
    Ffneu = [dat.F.FcellNeu{jj}(cellindex,:)];
    
    Coef='1/Regress';
    Ftmp = subtract_neurop(Ff,Ffneu,Coef);
    
    medF=medfilt1(Ftmp,dat.ops.imageRate*300,[],2,'truncate');
    Ftmp=Ftmp-medF;
    Ftmp=detrend(Ftmp','linear');
    
    [Ttmp,Ncell]=size(Ftmp);
    T(jj)=Ttmp;
    dF = cat(1,dF,Ftmp);
    
    
end
tmp = zeros(sum(T),Ncell);
 
percent_done = 0;
print_increment = 5;
fprintf(1, 'Estimating spikes from Ca trace ...\n');
parfor ii=1:length(cellindex)
%      fraction_done = 100 * (ii / length(cellindex));
%         if (fraction_done >= percent_done)
%             fprintf(1, '\b\b\b%d%%', percent_done);
%             percent_done = percent_done + print_increment;
%         end
%     fprintf('%d/%d done\n',jj,length(cellindex));
    dFtmp = dF(:,ii);
    dFtmp = dFtmp-min(dFtmp)+1;
    tmp(:,ii)=deconv_kh_002(dFtmp);
end
 fprintf(1, 'Done\n');
 StartInd = [1,cumsum(T)+1];StartInd = StartInd(1:end-1);
 EndInd   = [cumsum(T)];
 for jj=1:NumSubDirs
     index = StartInd(jj):EndInd(jj);
     dat.F.Spk{jj} = tmp(index,:)';
 end
fprintf('Saving %s...\n',FullProcSpkFile);

save(FullProcSpkFile,'-struct','dat');
