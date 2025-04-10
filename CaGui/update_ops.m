function ops = update_ops(ops_new,ops_old,varargin)
% usage: ops=update_ops(ops_new,ops_old);
% == input ==
% ops_new: newer options
% ops_old: oldrer options
% 
% if there is a same field in ops_old, ops_new, value in ops_new is used by default. 
% usage: ops=update_ops(ops_left,ops_right,'use_left');  % if there are same parameters, use left. (Default). 
% usage: ops=update_ops(ops_left,ops_right,'use_right'); % if there are same parameters, use right.  If the value on left is not on right, not used.
% ex)  
% ops_left.var1=1;
% ops_left.var2=2;
% ops_right.var1=10;
% ops_right.var3=3;
% ops=update_ops(ops_left,ops_right,'use_left')
% ops=update_ops(ops_left,ops_right,'use_right')
% 
% ex) new option: 'use_left_if_property_exist_in_right'
% ops is defined by default opt_right, but updated if the same property exists in left. 
% usage: ops=update_ops(ops_left,ops_right,'use_left_if_property_exist_in_right');
% 


if nargin>=3
    opt=varargin{1};
else
    opt='use_left';
end

FN_old = fieldnames(ops_old);
FN_new = fieldnames(ops_new);
switch opt
    case 'use_left'
        FN_only_in_old = setdiff(FN_old,FN_new);
        ops = ops_new;
        for ii=1:length(FN_only_in_old)
            ops.(FN_only_in_old{ii})=ops_old.(FN_only_in_old{ii});
        end
    case 'use_right'
        FN_only_in_new = setdiff(FN_new,FN_old);
        ops = ops_old;
        for ii=1:length(FN_only_in_new)
            ops.(FN_only_in_new{ii})=ops_new.(FN_only_in_new{ii});
        end
    case 'use_left_if_property_exist_in_right'
        FN_in_common = intersect(FN_new,FN_old);
        ops = ops_old;
        for ii=1:length(FN_in_common)
            ops.(FN_in_common{ii})=ops_new.(FN_in_common{ii});
        end
    otherwise
        error('Unknown option %s!',opt);
end
