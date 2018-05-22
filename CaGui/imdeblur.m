function [ out ] = imdeblur(img,h)
% imdeblur(img,h)
% img: gray scale image
% h:    filter image.
% ex)
% img = imread(file);
% h = fspecial('gaussian', [15 15], 3);
% 

img=im2double(img);

[m,n] = size(img);
[mb,nb] = size(h);

% output size 
mm = m + mb - 1;
nn = n + nb - 1;

out = ifft2(fft2(img,mm,nn)./(fft2(h,mm,nn)+eps));

out = out(1:m,1:n);