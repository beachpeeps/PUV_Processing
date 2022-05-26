function [hsig,hig,eta,CI_hsig, CI_hig,higlo,highi] = get_hsig_from_watlev(wlev,dt,nfft)
% [hsig,hig,eta,CI_hsig, CI_hig,higlo,highi] = get_hsig_from_wlev(wlev,dt,nfft)
%
% Inputs:   wlev: water level timeseries [nx,nt] matrix
%           dt: sampling time in seconds
%           nfft: size of the windows for spectral averaging
% Outputs:  hsig: Significant wave height at SS frequencies
%           hig: Significant wave height at IG frequencies
%           eta: mean water level
%           CI_hsig/hig: confidence intervals (at 95%)
%           higlo/highi: Sig wave height atlow and high IG frequency 
%
% Dependencies: 
%       get_spectrum
%           calculate_fft2  
%       getSWHebounds
% 
% TODO: This could probably use some edits for:
%   -   added flexibility to the frequency limits, 
%   -   confidence intervals
%   -   the way the window size is called?
%   -   option to do a subset of locations?


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

h = waitbar(0,'Calculating Hsig');

for i=1:nx
    w = waitbar(i/nx,h);
    x = wlev(i,:);
    
    [f, S, Spplo, Sppup, nens, dof] = get_spectrum(detrend(x)', nfft, 1/dt,0.05);
    nINC = find(f>=0.04 & f<=0.25);
    nIG = find(f>=0.004 & f<=0.04);
    nIGlo = find(f>0.004 & f<0.02);
    nIGhi = find(f>0.02 & f<0.04);
    df = f(2)-f(1);
    
    [ hsig(i), lb, ub, edof(i) ] = getSWHebounds( S(nINC), dof, 0.95, df);
    [ hig(i), lb_ig, ub_ig, edof_ig(i) ] = getSWHebounds( S(nIG), dof, 0.95, df);
    [ higlo(i), lb_iglo, ub_iglo, edof_ig(i) ] = getSWHebounds( S(nIGlo), dof, 0.95, df);
    [ highi(i), lb_ighi, ub_ighi, edof_ig(i) ] = getSWHebounds( S(nIGhi), dof, 0.95, df);
    
    CI_hsig(i,:) = [lb ub];
    CI_hig(i,:) = [lb_ig ub_ig];
    
    eta(i) = nanmean(x); % tricky to define at beach/water interface!!
    
end
close(h)