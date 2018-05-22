function ops = update_ops(ops_new,ops_old)
% usage: ops=update_ops(ops_new,ops_old);
% == input ==
% ops_new: newer options
% ops_old: oldrer options
% 
% if there is a same field in ops_old, ops_new, value in ops_new is used.

FN_old = fieldnames(ops_old);
FN_new = fieldnames(ops_new);

FN_only_in_old = setdiff(FN_old,FN_new);
ops = ops_new;
for ii=1:length(FN_only_in_old)
    ops.(FN_only_in_old{ii})=ops_old.(FN_only_in_old{ii});
end
