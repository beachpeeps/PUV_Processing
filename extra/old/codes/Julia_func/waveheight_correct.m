function [Hs,Hig,Tp] = waveheight_correct(p,dt,h,Tmin,Tmax)

% This function computes significant wave height from bottom pressure
% assuming the sensor is at depth h (meters).

% june 2007 Mark Merrifield? 

np = length(p);
np2 = floor(np/2);
freq = (1:np2)/(np*dt);

%Hs
[hk,kk]=dispr(h,freq);
%Kpp = 1./cosh(hk);
Kpp = cosh(hk);
k = find(freq >= 1/Tmin | freq <= 1/Tmax);
Kpp(k) = 0;

fp = fft(detrend(p))/np;
fp = fp(2:np2+1).*Kpp;
fp = 2*abs(fp).^2;
Hs = 4*sqrt(sum(fp));
[m,im] = max(fp);
Tp = 1./freq(im);

%Hig

Kpp = cosh(hk);
k = find(freq > 1/Tmax);
Kpp(k) = 0;

fp = fft(detrend(p))/np;
fp = fp(2:np2+1).*Kpp;
fp = 2*abs(fp).^2;
Hig = 4*sqrt(sum(fp));