function [AS, F, xfreq] = spectra(X, sfreq, nfft, win)

% [AS F xfreq] = spectra(X, sfreq, nfft, win)
%
% Obligatory inputs:
%	X     - signal matrix (fft is done along the first dimension)
%	sfreq - sampling frequency
%
% Optional inputs (defaults):
%	nfft = size(X,1); % length of samples in FFT (can be used for zero-padding)
%	win  = @hann; % @blackmanharris, @rectwin, @hann, @hamming, @blackman % window function applied to the data
%
% Outputs:
%	S     - power spectrum
%	AS    - amplitude spectrum
%	P     - phase spectrum
%	xfreq - x-axis for power/amplitude spectrum
%	fs    - step size of the frequency axis
%
% Description: The program computes the power, amplitude, and phase
% spectrum along the first dimension of a given matrix X. Gives only one half of the spectrum.
%
% References:
% van Drongelen W (2007) Signal processing for neuroscientists:
%	Introduction to the analysis of physiological signals. Amsterdam and
%	others: Academic Press.
%
% For infos on windowing see also: http://zone.ni.com/devzone/cda/tut/p/id/4844
% -----------------------------------------------------------------------
% B. Herrmann, M. Henry, Email: bherrmann@cbs.mpg.de, 2012-02-02

% check inputs
nInputs = nargin;
S = []; AS = []; P = []; xfreq = []; fs = [];
if nInputs < 2, fprintf('Error: Wrong number of inputs!\n'); return; end;
if nInputs < 3 || isempty(nfft), nfft = size(X,1); end;
if nInputs < 4 || isempty(win),  win = @hann; end;
if nfft < size(X,1); fprintf('Info: nfft < size(X,1)! Using nfft = size(X,1).\n'); nfft = size(X,1); end

% reshape multidimensional data into 2-d matrix
siz = size(X);
n = siz(1); nrest = siz(2:end);
X = reshape(X,[n prod(nrest)]);

% get window function and multiply it with the data
xwin = window(win,n);
X = bsxfun(@times,X,xwin);

% % compute FFT, F - complex numbers, fourier matrix
% F = fft(X,nfft,1);
% 
% compute cronux spectra
F = fft(X,nfft,1);


% compute power, amplitude, and phase spectra
AS = 2/sum(xwin) * sqrt(F.*conj(F)); % amplitude spectrum = abs(F)*2/sum(xwin) % sum(xwin) == n for rectangular window
S  = F.*conj(F)/sum(xwin); % power spectrum = (abs(F).^2)/sum(xwin); % sum(xwin) == n for rectangular window
P  = atan(imag(F)./real(F)); % phase spectrum

% get frequency axis
xfreq = (0:nfft/2)*(sfreq/nfft);
fs    = sfreq/nfft;

% shorten matrix by half
AS = AS(1:length(xfreq),:);
S  = S(1:length(xfreq),:);
P  = P(1:length(xfreq),:);
F  = F(1:length(xfreq),:);

% put data back into multi-d matrix
AS = reshape(AS,[size(AS,1) nrest]);
S  = reshape(S,[size(S,1) nrest]);
P  = reshape(P,[size(P,1) nrest]);
F  = reshape(F,[size(F,1) nrest]);

return