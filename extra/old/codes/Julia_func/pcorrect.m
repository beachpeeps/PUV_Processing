function P = pcorrect(pressure,Hz,fcutoff,burial)
% code takes in pressure at depth, measured with pressure sensor, and
% returns linearly corrected surface pressure time series.
test
if nargin<2
    fcutoff = 0.33;
end

% get mean wave height
H = nanmean(pressure);

% give error for psensors not under water table
if H<0
    error('Warning: negative H detected - sensor is likely over watertable: skipping hour');
end

%% Apply depth correction to frequencies below fcutoff
% % Take into frequency space

% infill nans, although nans should not exist...
pressure(isnan(pressure)) = nanmean(pressure);

%****** prewhiten pressure******
dp = diff(pressure);
pre = nanmean(pressure(1:100))*ones(600,1)+nanmean(dp(1:100)).*(1:600)'-nanmean(dp(1:100)).*600;
post = nanmean(pressure(end-100:end))*ones(600,1)-nanmean(dp(1:100)).*(1:600)';
ppad = [pre; pressure; post(1:end)];
%*******************************

%****create frequency vector****
nn=length(ppad);
ny = floor(nn/2);
f = Hz/2*linspace(-1,1,nn);

fY = fft(ppad);
fY = fftshift(fY);
%********************************

% locate frequencies subject to correction
cutoffind = find(abs(f) <= fcutoff & abs(f)>abs(f(ny)));

% do correction in frequency/wavenumber space
if H>0 % can only solve for k if water depth is greater than 0   
    ff = cutoffind;
    k(ff) = getk(f(ff),H);
    fY(ff) = fY(ff)'.*exp(k(ff).*burial).*cosh(k(ff)*H);    
end
%******************

%***bring back to time domain and remove padding
Psspad = real(ifft(ifftshift(fY)));
P = Psspad(601:end-600);
% need to solve the ringing problem at jumps/bores!!

