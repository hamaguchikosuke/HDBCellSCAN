mimg = handles.dat.ops.mimg1;
img=mimg-my_conv2(mimg,2,[1,2]);
nimg=img./my_conv2(sqrt(img.^2),4,[1,2]);
%%

figh=myfigure('ROI check');cla;
ax(1)=subplot(2,2,1);
imagesc(nimg);
% colormap gray;

celldiameter = 11;
cell_membrane = 2;
cell_edge  = 2;

hsize = celldiameter+cell_edge+cell_membrane;
h_boundary = zeros(hsize,hsize); 
h_boundary_ind = 1:hsize;
h_boundary(h_boundary_ind,h_boundary_ind)=fspecial('disk',(length(h_boundary_ind)-1)/2);

h_nuclei       = zeros(hsize,hsize);
h_nuclei_ind = (1+cell_membrane+cell_edge):(hsize-cell_membrane-cell_edge);
h_nuclei(h_nuclei_ind,h_nuclei_ind)=fspecial('disk',(length(h_nuclei_ind)-1)/2);

h_mem = zeros(hsize,hsize); 
h_mem_ind = (1+cell_membrane):(hsize-cell_membrane);
h_mem(h_mem_ind,h_mem_ind)=fspecial('disk',(length(h_mem_ind)-1)/2);


h = 2*h_mem-h_nuclei-h_boundary;
conved_nimg= conv2(nimg,h,'same');
ax(2)=subplot(2,2,2);
imagesc(conved_nimg);


%

nimg_rgb = repmat(mat2gray(nimg),[1,1,3]);
cell_center = zeros(size(nimg_rgb));
cell_center(:,:,1)=conved_nimg>0.8;
cell_center(:,:,2)=conved_nimg<-0.9;
ax(3)=subplot(2,2,3);
imagesc(tanh(nimg_rgb+cell_center));


% deconv_cell_center=zeros(size(nimg_rgb));
% deconv_cell_center(:,:,1)=imdeblur(double(cell_center(:,:,1)),fspecial('gaussian',[5 5],1));
% ax(4)=subplot(2,2,4);
% imagesc((1+tanh(nimg_rgb+deconv_cell_center))/2);
 
linkaxes(ax,'xy');