function [etaAa, etaAss, etaAig]= timeseries_ss_ig(eta, time)
%% Syntax

%% 
%% Author
% Athina Lange, SIO July 2019
%% 
[m,n] = size(eta);
timedat = (1:m)';
%% Frequency and wavenumber
if isdatetime(time) == 0 % if not datetime
   dt = nanmean(median(diff(time)),2) % in spc
else
    dt = seconds(nanmean(median(diff(time)),2)) % in spc
end
% nfft=size(eta,1);
% df = 2/(nfft-1);   % This is the frequency resolution
% nnyq = nfft/2 +1;
% frequency = [0:df:(nnyq-1)*df -(nnyq-2)*df:df:-df];

nfft = size(eta,1);
if isodd(nfft)
    nfft=nfft-1;
end

N = nfft;
fs = 1/dt;
frequency = [0:fs/N:fs/2 -fs/2+fs/N:fs/N:-fs/N];


%% Pressure response factor (inverse)
flim = [0.004 0.04 0.25]; % limits of IG and SS bands

%% Calculating eta
pdF = fft(eta,nfft);

% All
etaF = pdF; etaF(abs(frequency) > flim(3)) = 0; etaF(abs(frequency) < flim(1)) = 0;
etaA = ifft(etaF);
etaAa = etaA;

%SS
etaFss = pdF; etaFss(abs(frequency) > flim(3)) = 0; etaFss(abs(frequency) <= flim(2)) = 0;
etaAss = ifft(etaFss);

%IG
etaFig = pdF; etaFig(abs(frequency) > flim(2)) = 0; etaFig(abs(frequency) <= flim(1)) = 0;
etaAig = ifft(etaFig);

end