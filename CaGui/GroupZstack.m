function F = GroupZstack(F,n)
% Grouped Z-stack of F
% Usage: 
% n = 10; % average of n slices
% groupF = GroupZStack(F,n)
% 
% by KH 20170518
[H,W,T]=size(F);

if rem(T,n)~=0
    error('number of group %d must be a divisor of number of zstack %d',...
        n,T);
end


F = reshape(F,[H,W,n,T/n]);
F = squeeze(mean(F,3));