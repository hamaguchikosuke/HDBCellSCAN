function [ops, stat, res] = get_HDBCellScan(ops, U, Sv,clustrules)
% wrapper function of HDBCellScan to return res, stat, and ops
%
%

MinClust=clustrules.MinClust;

%L: labels
%prob: probability of being in the cluster
%M: max correlation to neighbor pixels
%S: neurpil 

[L,probabilities,M,S]  = HDBCellScan_20171102(U,Sv,MinClust);
% [L,probabilities,M,S]  = HDBCellScan_20180218(U,Sv,MinClust,S);
[Ly,Lx]=size(L);

stat    = get_stat_from_iclust(L,M);

res.iclust  = L;
res.M       = M;
res.S       = S;
res.lambda  = probabilities;
res.probabilities = probabilities;
%
res.Ly  = Ly;
res.Lx  = Lx;