function [f,fN,fstop,Sxxhat, b2hat, b2hat_complex,bphase, N,M,varx, Bhat] = bispectra_juliatrial(x,Fs,fmax,N)
% [f,fN,fstop,Sxxhat, b2hat, N,M,varx] = mybispectra(x,Fs,N,fmax)
% f = frequencies
% fN = Nyquist frequency
% fstop = index of largest freq to plot
% Sxxhat = variance preserving PSD
% b2hat = bicoherenceˆ2
% N = ensemble record length
% M = number of ensembles
% varx = variance of truncated timeseries



L = length(x)
fN = Fs/2;
if nargin<4
N = 128; % default ensemble record length
end

M = floor(L/N); % number of ensembles
%%
L = M*N % define new length

% split x up into individual records
xi = reshape(x(1:N*M),N,M);

% subtract mean from each record
xi = xi-ones(N,1)*mean(xi); % i don't like this.

%apply window to each record
w = hann(N);

for m=1:M
xw(:,m) = xi(:,m).*w;
end

% normalize the data after windowing it to the variance of original
% timeseries
varx = var(xi(:));
xw = xw/sqrt(var(xw(:))/varx);

% compute the fourier amplitudes
X = fft(xw);
% X = X(2:end,:); %ignore the DC component (why?? but okay...)


% Define frequency axis
df = Fs/N;
% f = Fs/2*linspace(0,1,N/2+1);
f = df:df:Fs/2;

% sum over all segments

% Energy Density, units: mˆ2/Hz (Avg Energy per frequency band)
Sxxhat = sum(abs(X.^2),2)/M;%/N;

% estimate the bispectrum
Bhat = nan(N/4,N/2);
de1 = nan(N/4,N/2);
de2 = nan(N/4,N/2);
bhat = nan(N/4,N/2);
bhat_complex = nan(N/4,N/2);

fstop = find(f<fmax,1,'last');
if mod(fstop,2)==1
    fstop = fstop-1;
end

% 5. estimate the bispectrum
Bhat = nan(fstop/2,fstop);
denom = nan(fstop/2,fstop);
b = nan(fstop/2,fstop);
b_complex = nan(fstop/2,fstop);
bphase = nan(fstop/2,fstop);

for l=1:fstop
    for k=1:floor(fstop/2)
        if k+l <= fstop && k <= l
            Bhat(k,l) = 1/L*sum(X(k,:).*X(l,:).*X(k+l,:),2); %skewness per frequency
            denom(k,l) = Sxxhat(k)*Sxxhat(l)*Sxxhat(k+l); %variance of skewness
            b(k,l) = abs(Bhat(k,l))/sqrt(denom(k,l)); %normalized bispectrum
            b_complex(k,l) = Bhat(k,l)/sqrt(denom(k,l));
            bphase(k,l) = atan2(imag(Bhat(k,l)),real(Bhat(k,l)));%biphase
        end
    end
end
b2hat = b/M; %average over all ensembles
b2hat_complex = b_complex/M; %average over all ensembles
Sxxhat = Sxxhat(1:N/2)*df; % variance preserving

