function Fout = subtract_neurop(Fcell,Fcellneu,varargin)
% subtract neuropil signal 
% F = Fcell-Coef*b(2)*Fcellneu
% where b(2) is the regression slope.
% Coef is empirically employed coefficient for subtraction, default is 0.5.
% 
% -- input ---
% Fcell: N x T matrix. Each row is fluorescence of ROI.
% Fcellneu: N x T matrix. Each row is fluorescence of neuropil mask.
% Coef: can be empty. 
% 
% -- output --
% F = Fcell-Coef*b(2)*Fneurop
% where b(2) is the regression slope.
% 
% 
% ex) 
% T = 1000;
% base = randn(1,T);
% TrueF = 0.2*randn(1,T);
% Fcell=     TrueF+base;
% Fcellneu = 10*base+randn(1,T);
% Fout = subtract_neurop(Fcell,Fcellneu);
% 
% plot(TrueF,Fout)
% 
% by KH 20171124

if nargin>=3
    Coef = varargin{1};
else
    Coef = 0.5;
end

[N,T]=size(Fcell);
[N1,T1]=size(Fcellneu);

if N~=N1 | T ~= T1
    error('size of Fcell and Fcellneu must be the same');
end

Fout = zeros(size(Fcell));
for nn=1:N
    y = Fcell(nn,:)';
    x = [ones(T,1), Fcellneu(nn,:)'];
    % b=regress(y_trace,x);
    b = (x'*x)\(x'*y);
    Coef = 0.6;
    
    Fout(nn,:)=y-Coef*b(2)*Fcellneu(nn,:)';
    minFout = min(Fout(nn,:));
    
    if minFout<0
        Fout(nn,:)=  Fout(nn,:)-minFout*2;
    end
end