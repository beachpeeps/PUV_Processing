function [Pss,posEta] = depthCorrect(PUV,tnum)

P = PUV(tnum).P(:);
U = PUV(tnum).U(:);
%%


%get height above sand of pressure sensor
% cd ~/Agate/mfiles/PUVprocessing
% [doffp,doffu,~] = get_doffs(8,filename);
% cd(currdir)

% load ~/Agate/mat/Pmean pres_barometric time
% ind = knnsearch(time', filename);
% baro = pres_barometric(ind);
% P = P-baro;

doffp = 0.7; %70 cm for deployment 2
doffu = doffp+.22;
Fs = 2; %sampling frequency in hertz

depth = mean(P) + doffp;  % This is the hourly averaged water depth
P = detrend(P);
U = detrend(U);

N = length(P);
if mod(N,2)==1
    P=P(1:end-1);
    U=U(1:end-1);
    N = length(P);
end

%%

% 
NN = floor(N/2);
f = Fs/2*linspace(-1,1,N);
%%
fY = fft(P);
fY = fftshift(fY);

fU = fft(U);
fU = fftshift(fU);

fcutoff = 0.25;
cutoffind = find(abs(f) <= fcutoff);
cutoffind = find(abs(f) <= fcutoff & abs(f)>abs(f(NN)));
%depth correction
    for ff = cutoffind
        %currdir = pwd;
        %cd ~/Agate/mfiles
        k = ksolve(depth,2*pi*f(ff));
        fY(ff) = fY(ff)*cosh(k*depth)./cosh(k*doffp);
        fU(ff) = fU(ff)*cosh(k*depth)./cosh(k*doffu);

        %cd(currdir)
    end

fY = ifftshift(fY);
Pss = real(ifft(fY));

omega = 2*pi*f; g= 9.81;
% cdir = pwd;
% cd ~/Agate/mfiles/PUVprocessing/
k = get_wavenumber( omega , depth);
% cd(cdir)
convert = omega./(g*k);  % record conversion coefs to check them
%Turning U velocities into pressure
fU = fU.*convert';

fU = ifftshift(fU);
Upres = real(ifft(fU));

%Only incoming waves, only works for IG Shallow Water
posEta = 0.5*(Pss+Upres);