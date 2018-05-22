 function out=myfigure(opt,varargin)
% make a figure panel in similar way that figure() works, 
% but it also accepts texts as an identifier.
% ex) myfigure('experiment1')
% ex) myfigure('experiment','-regexp'); % you can search figure which name
% starts from "experiment", and rename it as "experiment".

if nargin==2
    findobj_option=varargin{1};
else
    findobj_option=[];
end

if isscalar(opt)
    fig1=figure(opt);
elseif ischar(opt)
    
    if isempty(findobj_option)
        fig1=findobj('Name',opt);
    else
        fig1=findobj(findobj_option,'Name',opt);
    end
    
%     fig1
    
    if isempty(fig1)
        num=ceil(1000*rand(1))+1000;
        fig1=figure(num);
        set(fig1,'Name',opt,'NumberTitle','Off');
    else
        figure(fig1(1));
    end
end

out=fig1(1);