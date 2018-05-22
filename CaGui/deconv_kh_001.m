%% test deconvolution of ca image trace to spike trains
load('Fcell_example.mat');
%%
myfigure('Deconv test');clf;
F=double(Fcell{1}(231,1:1000))';
Fmin = min(F);
Fmax = max(F);
F = (F-Fmin)/(Fmax-Fmin);
h_data=plot(F);hold on;


Fs = 33;
SensorTau = 1 * Fs;
alp = 1;                                           % rescalouing factor. 
bet = median(F);                            % estimate of baseline
gam = (1-1/Fs);
sig  = median(abs(F-bet))/1.4785; % estimate of noise from median absolute deviation.
sig2 = sig^2;
lambda = 1;                                    % estimate of poisson firing rate.
T = length(F);
%%  init parameters
% make bidiagonal matrix to compute C_{t}-gam*C_{t-1}.
gaus = gausswin(15); gaus = gaus(:)/sum(gaus);
C = conv2(F,gaus,'same');
h_C = plot(C,'r');hold on;
n = M*F;
h_n = plot(n,'k');
index_i = [1:T-1, 1:T-1];
index_j = [1:T-1,2:T];
s           = [-gam * ones(1,T-1),ones(1,T-1)];
M = sparse(index_i,index_j,s,T-1,T);
M = 
lmd = lambda * ones(T-1,1);
H2 = speye(T);
d0 = find(H2); % diagonal index 

ss = 0.001; % step size
% solve argmax_MC -1/(2*sig^2) abs(F-C-bet).^2-MC'*lmd + z * log(MC'*ones(T,1));
% z =0.001;
 Mlmd =  (1-gam)*lambda * ones(T,1);
%%
for z=[0.01,0.001,0.0001,0.00001]
    for ii=1:1000
        MC = M*C;
        OneOverMC = 1./(MC);
        g = -(F-C-bet)/sig2 + ones(T,1)*lambda*(1-gam) - z*M'*OneOverMC; % gradient
%         g = -(F-C-bet)/sig2 + Mlmd - z*M'*OneOverMC; % gradient, continuous boundary condition.
        H2(d0) =  OneOverMC.^2;
        H = speye(T)/sig2 - z * M' * H2 * M;
%         H = speye(T)/sig2 - z * M' *( M ./(repmat(OneOverMC.^2,1,T)));    % Hessian
        d=H\g;
        C =  C + ss * d;
        n = M*C;
        n(n<0)=0;
        lambda = T*Fs/sum(n);
        
%         Pzold =  -(F-C-bet).^2/(2*sig2)-MC'*lmd + z*sum(log(MC));
        
%         Pz = -(F-C-bet).^2/(2*sig2)-MC'*lmd + z*sum(log(MC));
        
        bet = max(0,mean(F-C));
%         bet = min(F);
        sig2 = mean((F-C-bet).^2);
%         lmd = lambda * ones(T-1,1);
        Mlmd =  (1-gam)*lambda * ones(T,1);
        set(h_C,'YData',C);
        set(h_n,'YData',n);
        if mod(ii,10)==0, drawnow; end
    end
end