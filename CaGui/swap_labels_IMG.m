function Lout = swap_labels_IMG(L,before_label,after_label)
%  L = swap_labels_IMG(L,original_label,swap_label);
%  change the label in the roi label matrix
%  
% -- input --
% L: HxM integer matrix. 0 is background. integers are ROI labels.
% original_label: 1xN or Nx1 vector, indicating the label before swap.
% swap_label    : 1xN or Nx1 vector, indicating the label after swap.
% 
% -- output --
% L: HxM integer matrix after swapping labels.
% ex)
% L=[0 0 1; 0 2 0; 3 0 0]
% before_label = [1:3];
% after_label = [11,21,31];
%  L = swap_labels_IMG(L,before_label,after_label)

N_before_label = length(before_label);
N_after_label = length(after_label);

if N_before_label ~= N_after_label
    error('before and after_label size must be the same');
end

Lout = L;
for ii=1:N_before_label
    Lout(find(L==before_label(ii)))=after_label(ii);
end