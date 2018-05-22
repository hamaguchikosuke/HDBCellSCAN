function [n,C,P]=deconv_kh_003(F,varargin)
% [n,C,P]=deconv_kh_003(F,dim)
% ==== Input ===
% F : T x N (default) or N x T matrix where T is the length of time, N is the number of cells. 
% dim: dimension of time. Default is 1. 
% 
% ==== Output ====
% n:  estimated spikes 
% C:  estimated Ca trace
% P:  estimated parameters 
% 
% This program calculates MAP estimate of spike train from fluorescence dynamics 
% based on Vogelstein et al.,  J. Neurophysiol., 2010. (Fast_oopsi)
% 
% ex) test deconvolution of ca image trace to spike trains
% 
% load('Fcell_example.mat');
% F=double(Fcell{1}(231,1:1000))';
% Fmin = min(F);
% Fmax = max(F);
% F = (F-Fmin)/(Fmax-Fmin);
% [n,C,P]=deconv_kh_002(F);
% plot(F,'c'); hold on; plot(P.alp*C+P.bet,'b'); plot(n,'k');
% 
%%

if nargin>=2
    dim = varargin{1};
else
    dim = 1;
end

switch class(F)
    case 'single'
        F = double(F);
    otherwise
end

if dim==1
   % do nothing.
elseif dim==2
   F = F'; % Now TxN matrix
end

[T,nCell]=size(F);
n=nan(size(F));
C = nan(size(F));
P = nan(nCell,3);
for ii=1:nCell
    if all(F(ii,:)==0) % not an roi
        % do nothing
    else
        [ntmp,Ctmp,Ptmp]=local_deconv(F(:,ii));
        n(:,ii)=ntmp;
        C(:,ii)=Ctmp;
        P(ii,:)=P;
    end
end
end
% F=double(Fcell{1}(231,1:1000))';
% F = detrend(F);
% F = F-min(F)+eps;
function [n,C,P]=local_deconv(F)
persistent cnt
if isempty(cnt)
    cnt=1;
else
    cnt = cnt+1;
end
do_plot =1;

Fs = 33;
SensorTau = 1 * Fs;
qF = quantile(F,[0.5,0.05]);                
alp = qF(1);                                        % rescalouing factor.
bet = qF(2);                                       % estimate of baseline
gam = (1-1/SensorTau);
sig  = median(abs(F-bet))/1.4785; % estimate of noise from median absolute deviation.
sig2 = sig^2;
lambda = 10;                                    % estimate of poisson firing rate.
T = length(F);
%%  init parameters
% make bidiagonal matrix to compute C_{t}-gam*C_{t-1}.
% gaus = gausswin(15); gaus = gaus(:)/sum(gaus);
% C = conv2(F,gaus,'same');
C = bet * ones(T,1);

index_i = [2:T, 1:T];
index_j = [1:T-1,1:T];
s           = [-gam * ones(1,T-1),1-gam, ones(1,T-1)]; % using C(0)=1/gam * C(1) boundary condition.
M = sparse(index_i,index_j,s,T,T);

n = M*F;

lmd = lambda * ones(T-1,1);
H2 = speye(T);
d0 = find(H2); % diagonal index

s = 1; % step size
% solve argmax_MC -1/(2*sig^2) abs(F-C-bet).^2-MC'*lmd + z * log(MC'*ones(T,1));

Mlmd =  (1-gam)*lambda * ones(T,1);

if (do_plot) 
    figure(400);clf;
    h_data=plot(F);hold on;
    h_C = plot(C,'r');hold on;
%     h_d = plot(NaN(T,1),'c');hold on;
%     h_g = plot(NaN(T,1),'g');hold on;
    h_n = plot(n,'k');
    title(sprintf('cell %d',cnt));
end

%%
for z= (0.1).^[0:13]
    
    ss=1;
    d =1;
    cnt = 1;
    while  norm(d)>0.05 && ss>1.e-3
        cnt = cnt+1;
        n = M*C;
        OneOverMC = 1./n;
        D = (F-alp*C-bet); % difference vector
        g = alp*D/sig2 - ones(T,1)*lambda*(1-gam) + z*M'*OneOverMC; % gradient
        %         g = -(F-C-bet)/sig2 + Mlmd - z*M'*OneOverMC; % gradient, continuous boundary condition.
        H2(d0) =  OneOverMC.^2;
        H =  alp^2*speye(T)/sig2 + z * M' * H2 * M;
        %         H = speye(T)/sig2 - z * M' *( M ./(repmat(OneOverMC.^2,1,T)));    % Hessian
        d=H\g;
        
        
        post = -1/(2*pi*sig2)*sum(D.^2) - n' * Mlmd + z * sum(log(n));
        post1 = post - 1;
        
        hit = -n./(M*d);                % step within constraint boundaries
        hit(hit<0)=[];                  % ignore negative hits
        
        if any(hit<1)
            ss = min(1,0.999*min(hit(hit>0)));
        else
            ss = 1;
        end
        
        while  post1 < post-1e-7
            C1 =  C + ss * d;
            D = F- alp*C1-bet;
            n = M*C1;
            post1 = -1/(2*pi*sig2)*sum(D.^2) - n' * Mlmd + z * sum(log(n));
            ss = ss/2;
            if ss<1.e-20, break; end
        end
        
        post = post1;
        C = C1;
        
        if (do_plot)
            set(h_C,'YData', alp*C+bet);
            set(h_n,'YData',n);
%             set(h_g,'YData',alp*CC+bet);
%             set(h_d,'YData',d);
%             pause(0.5);
        end
    end
%     fprintf('z=%2.0g,\tcnt=%d ss=%2.1g\n',z,cnt,ss)
    
         % estimate alp and bet from linear regression
        nsort = sort(n);
        nn=zeros(T,1);
%         nthr= 0.001;
        nthr = nsort(ceil(0.97*T));
        nn(n>nthr)=1;
        CC = filter(1,[1,-gam],nn);
        A = [CC, ones(T,1)];
        X = A\F;
        alp = X(1);
%         bet = X(2);
        bet = qF(2);
        sig2 = mean(D'*D)/T;
        
end


if nargout>=3
    P.alp = alp;
    P.bet = bet;
end
end