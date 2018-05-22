[F,info]=BigTiffReader;
FullFilename = info(1).Filename;
[FilePath,FileName,ext]=fileparts(FullFilename);
%% 
info(1).ImageDescription

[W,H,D]=size(F);
depth=3*[1:D];

myfigure('Movie');
ref=F(:,:,60);
imshow(ref);caxis([0 2^12]);

FullVideoName = fullfile(FilePath,[FileName,'.avi']);
v= VideoWriter(FullVideoName,'Uncompressed AVI');
%%
open(v);
for dd=1:D
%     tmp=imhistmatch(F(:,:,dd),ref);
tmp=F(:,:,dd);
    imshow(tmp);caxis([0 2^12]);
    depthtext = sprintf('%3dum',depth(dd));
    textH=text([max(xlim)-0.15*diff(xlim)],max(ylim)-0.05*diff(ylim),depthtext,'Color',[1 1 1]);
    drawnow;
    frame = getframe;
    writeVideo(v,frame);
end
close(v);
