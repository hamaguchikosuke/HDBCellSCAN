%  how to register two fluorescent plane 

%% first, read two big stack of mClover image
[IMG,info]=BigTiffReader;
mean_IMG_Fixed=mean(IMG,3);
IMG_filename_Fixed = info(1).Filename;

[IMG,info]=BigTiffReader;
mean_IMG_Target=mean(IMG,3);
IMG_filename_Target = info(1).Filename;

clear IMG;
%% prepare an analysis directory
cnt = 1;
while IMG_filename_Fixed(cnt)==IMG_filename_Target(cnt), cnt=cnt+1; end
[CommonPath,~,~] = fileparts(IMG_filename_Fixed(1:cnt-1));

AnalysisPath = fullfile(CommonPath,'Analysis');
if ~exist(AnalysisPath,'dir')
mkdir(AnalysisPath);
end
%% 
myfigure('Fixed');
imagesc(mean_IMG_Fixed);
title('Fixed');
colormap gray;
[~,TiffName_Fixed,ext] = fileparts(IMG_filename_Fixed);
printname = sprintf('mean_%s',TiffName_Fixed);
print(gcf,'-dtiff',fullfile(AnalysisPath,printname))


myfigure('Target');
imagesc(mean_IMG_Target);
title('Target');
colormap gray;

[~,TiffName_Target,ext] = fileparts(IMG_filename_Target);
printname = sprintf('mean_%s',TiffName_Target);
print(gcf,'-dtiff',fullfile(AnalysisPath,printname))

%% register 
[optimizer,metric] = imregconfig('multimodal');
optimizer.InitialRadius = optimizer.InitialRadius/3.5;
optimizer.MaximumIterations = 300;

moved_target = imregister(mean_IMG_Target,mean_IMG_Fixed, 'affine', optimizer, metric);

myfigure('Registered');clf;
imshowpair(moved_target, mean_IMG_Fixed,'falsecolor','ColorChannels','red-cyan');
title('Cyan: fixed, Red: moved target');

filename = [TiffName_Fixed, '_', TiffName_Target];
print(gcf,'-dtiff', fullfile(AnalysisPath,[filename,'.tiff']) );
