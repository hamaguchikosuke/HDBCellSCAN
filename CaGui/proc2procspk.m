function proc2procspk(FullProcFile)
%% load proc and add spike convolution data
%  proc2procspk(FullProcFile);
%  FullProcFile: full path and filename of *_proc.mat
%  Output is _procSpk.mat
% 
%  by KH 20171121


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

    for jj=1:NumSubDirs
        
        cellindex=find(dat.cl.selected);
         if max(cellindex)>size(dat.F.Fcell{jj},1)
             error('it seems that proc file is not finalized (cellindex size mismatched with dat.F.Fcell size)');
         end
        Ff=[dat.F.Fcell{jj}(cellindex,:)];
        Ffneu = [dat.F.FcellNeu{jj}(cellindex,:)];
        
        F = subtract_neurop(Ff,Ffneu);
        
        [Ncell,T]=size(F);
        % x = 0;
        % waitH=waitbar(x,'Deconvolving Ca into spikes');
        tmp = zeros(T,Ncell);
        parfor ii=1:length(cellindex)
%             ii
            tmp(:,ii)=deconv_kh_002(F(ii,:)');
        end
        fprintf('%d/%d done\n',jj,NumSubDirs);
        dat.F.Spk{jj} = tmp';
    end
    fprintf('Saving %s...\n',FullProcSpkFile);
  
    save(FullProcSpkFile,'-struct','dat');
  