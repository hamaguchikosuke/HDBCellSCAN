function [skewF, SNratio,C_of_zF, C_of_dzF]=get_F_dependent_stats(F)

skewF = skewness(F,0,2);

zF = zscore(F,0,2);
d = 10;
dzF = zscore(zF(:,(1+d):end)-zF(:,1:end-d),0,2);

% cellid = find(h.dat.cl.selected); % to visualize correlation, 
% h.dat.cl.C_of_zF=zF(cellind,:)*zF(cellind,:)'/size(zF,2);
% would be easier to see.  
C_of_zF=zF*zF'/size(zF,2);
C_of_dzF=dzF*dzF'/size(dzF,2);


%% Here is a new code to add signal/noise ratio by calculating the power ratio below 1Hz and above. 
FPS= 30;
Tlen = 100*FPS;
% Multitaper Time-Frequency Power-Spectrum (power spectrogram)
% function A=mtpsg(x,nFFT,Fs,WinLength,nOverlap,NW,nTapers)
% x : input time series
nFFT = 2^nextpow2(Tlen); %number of points of FFT to calculate (default 1024)
Fs = FPS; %sampling frequency (default 2)
WinLength = nFFT; %length of moving window (default is nFFT)
nOverlap = 0;%nFFT/2; %overlap between successive windows (default is WinLength/2)
NW = 3; %time bandwidth parameter (e.g. 3   or 4), default 3
nTapers = 2*NW-1; % nTapers = number of data tapers kept, default 2*NW -1
Detrend= 1; 

SNratio = nan(1,size(F,1));

waitH = waitbar(0,'Computing power spectrum ...');
for ii=1:size(zF,1)
[yo, fo]=my_mtcsg(zF(ii,:),nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers);
myo=mean(abs(yo),2);
% plot(fo,myo);

ind=fo<1; 
sig1Hz_power=sum(myo(ind));
noise_power=sum(myo(~ind));
SNratio(ii)=sig1Hz_power/noise_power;
waitbar(ii/size(zF,1),waitH);
end
close(waitH);
SNratio(isnan(SNratio))=1;
SNratio = log(SNratio(:));