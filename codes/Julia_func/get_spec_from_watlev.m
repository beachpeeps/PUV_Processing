function [Spec] = get_spec_from_watlev(wlev,dt,nfft)
% [Spec] = get_spec_from_wlev(wlev,dt,nfft)
%
% Inputs:   wlev: water level timeseries [nx,nt] matrix
%           dt: sampling time in seconds
%           nfft: size of the windows for spectral averaging
% Outputs:  Spec structure with f, S, and confidence limits
%
% Dependencies: 
%       get_spectrum
%           calculate_fft2  
% 
% TODO: update documentation


% copyright 2019 Julia W Fiedler jfiedler@ucsd.edu

% read in water level, it should be nx x nt
[nx,nt] = size(wlev);
if nx > nt
    disp('wlev variable should be (distance,time): flipping now')
    wlev = wlev';
end
%TODO: add options to do a subset of x locations

hsig = nan(nx,1);
hig = nan(nx,1);
higlo = nan(nx,1);
highi = nan(nx,1);
CI_hsig = nan(nx,2);
CI_hig = nan(nx,2);
eta = nan(nx,1);


for i=1:nx
    x = wlev(i,:);
    
    [f, S, Spplo, Sppup, nens, dof] = get_spectrum(detrend(x)', nfft, 1/dt,0.05);

    Spec(i).f = f;
    Spec(i).S = S;
    Spec(i).Spplo = Spplo;
    Spec(i).Sppup = Sppup;
    Spec(i).nens = nens;
    Spec(i).dof = dof;
    
    
end
