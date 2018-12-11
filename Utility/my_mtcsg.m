function [yo, fo, to, varargout]=my_mtcsg(varargin)
% Multitaper Time-Frequency Cross-Spectrum (cross spectrogram)
% function [yo, fo, to]=mtcsg(x,nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers);
% function [yo, fo, to, stat]=mtcsg(x,nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers);
% x : input time series
% nFFT = number of points of FFT to calculate (default 1024)
% Fs = sampling frequency (default 2)
% WinLength = length of moving window (default is nFFT)
% nOverlap = overlap between successive windows (default is WinLength/2)
% NW = time bandwidth parameter (e.g. 3 or 4), default 3
% nTapers = number of data tapers kept, default 2*NW -1
%
% output yo is yo(f, t)
%
% If x is a multicolumn matrix, each column will be treated as a time
% series and you'll get a matrix of cross-spectra out yo(f, t, Ch1, Ch2)
% NB they are cross-spectra not coherences. If you want coherences use, mtcohere
% 
% stat.varX     = variance of X within the moving window, nFFTChunks x nChannels matrix 
% stat.meanX    = mean of X within the moving window, nFFTChunks x nChannels matrix


% Original code by Partha Mitra - modified by Ken Harris
% Also containing elements from specgram.m

% Kosuke Hamaguchi modified original mtcsg from KenToolbox from 
% http://qneuro.rutgers.edu/Software.html

% default arguments and that
[x,nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers,nChannels,nSamples,nFFTChunks,winstep,select,nFreqBins,f,t] = my_mtparam(varargin);

% allocate memory now to avoid nasty surprises later
y=complex(zeros(nFreqBins, nFFTChunks, nChannels, nChannels)); % output array
Periodogram = complex(zeros(nFreqBins, nTapers, nChannels, nFFTChunks)); % intermediate FFTs
Temp1 = complex(zeros(nFFT, nTapers, nFFTChunks));
Temp2 = complex(zeros(nFFT, nTapers, nFFTChunks));
Temp3 = complex(zeros(nFFT, nTapers, nFFTChunks));
eJ = complex(zeros(nFFT, nFFTChunks));

% calculate Slepian sequences.  Tapers is a matrix of size [WinLength, nTapers]
[Tapers V]=dpss(WinLength,NW,nTapers, 'calc');

% New super duper vectorized alogirthm
% compute tapered periodogram with FFT 
% This involves lots of wrangling with multidimensional arrays.

if nargout==4
    stat.varX   = zeros(nFFTChunks,nChannels);
    stat.meanX  = zeros(nFFTChunks,nChannels);
end
    
TaperingArray = repmat(Tapers, [1 1 nChannels]);
for j=1:nFFTChunks
	Segment = x((j-1)*winstep+[1:WinLength], :);
	if (~isempty(Detrend))
		Segment = detrend(Segment, Detrend);
	end;
    
	SegmentsArray = permute(repmat(Segment, [1 1 nTapers]), [1 3 2]);
	TaperedSegments = TaperingArray .* SegmentsArray;
						
	fftOut = fft(TaperedSegments,nFFT);
	Periodogram(:,:,:,j) = fftOut(select,:,:); %fft(TaperedSegments,nFFT);
  
    if nargout==4
        stat.varX(j,:) = var(Segment,1);
        stat.meanX(j,:) = mean(Segment,1);
    end
end	
	
% Now make cross-products of them to fill cross-spectrum matrix
for Ch1 = 1:nChannels
	for Ch2 = Ch1:nChannels % don't compute cross-spectra twice
		Temp1 = squeeze(Periodogram(:,:,Ch1,:));
		Temp2 = conj(squeeze(Periodogram(:,:,Ch2,:))); 
		Temp3 = Temp1 .* Temp2;
		eJ=sum(Temp3, 2);
%         size(eJ)
		y(:,:, Ch1, Ch2)= eJ/nTapers;
		
		% for off-diagonal elements copy into bottom half of matrix
		if (Ch1 ~= Ch2)
			y(:,:, Ch2, Ch1) = conj(eJ) / nTapers;
        end     
	end
end



% we've now done the computation.  the rest of this code is stolen from
% specgram and just deals with the output stage

if nargout == 0
	% take abs, and use image to display results
    newplot;
    for Ch1=1:nChannels, for Ch2 = 1:nChannels
    	subplot(nChannels, nChannels, Ch1 + (Ch2-1)*nChannels);
	    if length(t)==1
	        imagesc([0 1/f(2)],f,20*log10(abs(y(:,:,Ch1,Ch2))+eps));axis xy; colormap(jet)
	    else
	        imagesc(t,f,20*log10(abs(y(:,:,Ch1,Ch2))+eps));axis xy; colormap(jet)
	    end
	end; end;
    xlabel('Time')
    ylabel('Frequency')
elseif nargout == 1,
    yo = y;
elseif nargout == 2,
    yo = y;
    fo = f;
elseif nargout == 3,
    yo = y;
    fo = f;
    to = t;
elseif nargout == 4,
    yo = y;
    fo = f;
    to = t;
    varargout{1}=stat;
end

% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu