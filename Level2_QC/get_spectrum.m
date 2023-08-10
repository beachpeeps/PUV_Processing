function [fm, Spp, Spplo, Sppup, nens, dof] = get_spectrum(P, nfft, Hz,alpha)
% [fm, Spp, Spplo, Sppup, nens, dof] = get_spectrum(P, nfft, fs,alpha)
% Makes a spectrum for a single point only
%
% Input: 
%   P: timeseries to be analyzed (should be size (nt,1))
%   nfft: window size for spectral averaging
%   Hz: sampling rate (Hz)
%   alpha: for confidence intervals (alpha = 0.05 is standard 95% Conf Int)
%
% Output:
%   fm: frequency
%   Spp: Spectra (inputUnits^2/Hz)
%   Spplo/Sppup: lower and upper bounds of confidence intervals
%   nens: number of ensembles used
%   dof: degrees of freedom

% copyright 2019 J Fiedler jfiedler@ucsd.edu


[m,n] = size(P);
if m<n
    P = P';
end


[Amp,nens] = calculate_fft2(P(:,1),nfft,Hz);
nens2 = 0;

% TODO: FIX THIS, FOR NOW only use one-dim P
% if n==2
% [Amp2,nens2] = calculate_fft2(P(:,2),nfft,fs);
% Amp = [Amp; Amp2];
% end

df = Hz/(nfft-1);   % This is the frequency resolution
nnyq = nfft/2 +1;
fm = [0:nnyq-1]*df;


Spp = mean( Amp .* conj(Amp) ) / (nnyq * df);  % need to normalize by freq resolution
Spp = Spp(1:nnyq);

Spp=Spp(:);

nens = nens+nens2;
% Confidence Intervals (CI)
dof = 8/3*nens; % for Hanning window

a = dof.*sum(Spp).^2./(sum(Spp.^2)); % effective dof
adof = a/(1+2/dof); % unbiased edof

chi2 = [chi2inv(alpha/2,dof) chi2inv(1-alpha/2,dof)];
CI_chi2 = [(1/chi2(1)).*(dof*Spp) (1/chi2(2)).*(dof*Spp)];
Spplo = CI_chi2(:,1);
Sppup = CI_chi2(:,2);


end
