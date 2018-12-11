% helper function to do argument defaults etc for mt functions
% my_mtparam(x,nFFT,sF,WinLength,nOverlap,NW,Detrend,nTapers);
% 
% x is the input matrix, each column is the signal.
% <Default parameters>; 
% nFFT = 1024;  % fft 
% Fs = 2;       % sampling rate
% WinLength = nFFT; 
% nOverlap = WinLength/2; 
% NW = 3;
% Detrend = ''; 
% nTapers = 2*NW -1; 
% 
%

% KH modified mtparam from KenToolbox.
% 20090916

function [x,nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers,nChannels,nSamples,nFFTChunks,winstep,select,nFreqBins,f,t] = mtparam(P)

nargs = length(P);

x = P{1};
if (nargs<2 | isempty(P{2})) nFFT = 1024; else nFFT = P{2}; end;
if (nargs<3 | isempty(P{3})) Fs = 2; else Fs = P{3}; end;
if (nargs<4 | isempty(P{4})) WinLength = nFFT; else WinLength = P{4}; end;
if (nargs<5 | isempty(P{5})) nOverlap = WinLength/2; else nOverlap = P{5}; end;
if (nargs<6 | isempty(P{6})) NW = 3; else NW = P{6}; end;
if (nargs<7 | isempty(P{7})) Detrend = ''; else Detrend = P{7}; end;
if (nargs<8 | isempty(P{8})) nTapers = 2*NW -1; else nTapers = P{8}; end;

% Now do some compuatations that are common to all spectrogram functions

winstep = WinLength - nOverlap;


nChannels = size(x, 2);
nSamples = size(x,1);

% check for column vector input
if nSamples == 1 
	x = x';
	nSamples = size(x,1);
	nChannels = 1;
end;

% calculate number of FFTChunks per channel
% nFFTChunks = round(((nSamples-WinLength)/winstep));
if nSamples==WinLength
    nFFTChunks = 1;
else
    nFFTChunks = ceil((nSamples-WinLength)/winstep);
end

% turn this into time, using the sample frequency
t = winstep*(0:(nFFTChunks-1))'/Fs;

% set up f and t arrays
if ~any(any(imag(x)))    % x purely real
	if rem(nFFT,2),    % nfft odd
		select = [1:(nFFT+1)/2];
	else
		select = [1:nFFT/2+1];
    end	
else
	select = 1:nFFT;
end
nFreqBins = length(select);
f = (select - 1)'*Fs/nFFT;
