function RGB=ROI_gem_img(labels,probabilities,M,varargin)
% generate beautiful ROI image.
% 
% usage1)
% RGB = ROI_gem_img(labels,probabilities,M);
% 
% -- input --
% labels:  Ly x Lx integer matrix. Each ROI has a distinct label number as
% an integer.
% probabilities: Ly x Lx matrix with [0 1] values. Likelihood of each pixel belonging to its cluster.
% This one is calculated from HDBScan. This determines the saturation (0 is gray, 1 is vivid color).
% M :     Ly x Lx matrix with [0 1] values. This is a normalized, mean image
% 
% -- output --
% RGB: Ly x Lx x 3 RGB image.
% 
% 
% usage2)
% RGB = ROI_gem_img(labels,probabilities,M,opt);
% opt: options
% opt.hue_range : range of hues, ex) [0.2 0.21] is light green, [0 0.01]
% is deep red, [0.7 0.71] is blue. 
% opt.highlight_labels: a vector of indices. ROI belongs to this indices are colored with opt.highlight_color. 
% opt.highlight_color : a vector of indices to color specific ROIs. 
% by Kosuke Hamaguchi 20171117

if nargin>=4
    opt = varargin{1};
else
    opt = [];
end

opt.hue_range        = getOr(opt, 'hue_range',[0.2 0.8]);
opt.highlight_labels = getOr(opt, 'highlight_labels', []);
opt.highlight_color  = getOr(opt, 'highlight_color',  []);

if isscalar(opt.highlight_color)
    opt.highlight_color = opt.highlight_color*ones(length(opt.highlight_labels),1);
end
L = labels;
[Ly,Lx]=size(labels);
L(L==0)=max(unique(labels))+1; % usually, zero is the background. 

rng('default');
mid_color = mean(opt.hue_range);
color_range = diff(opt.hue_range)/2;
rand_hue  = mid_color+color_range*2*(rand(1,max(unique(L)))-0.5);
rng('shuffle');
rand_hue= [rand_hue,0];
Hue = rand_hue(L);
Sat = reshape(probabilities,Ly,Lx);
Val = M;

for ii=1:length(opt.highlight_labels)
   Hue(labels==opt.highlight_labels(ii))=opt.highlight_color(ii); % change color
end
    
RGB=hsv2rgb(cat(3,Hue,Sat,Val));