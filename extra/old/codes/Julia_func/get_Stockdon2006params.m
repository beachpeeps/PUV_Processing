function [R2,Sig,Sinc,eta,sqrtHoLo,S,xi] = get_Stockdon2006params(Ho,Lo,beta)
% [R2,Sig,Sinc,eta,sqrtHoLo] = get_Stockdon2006params(Ho,Lo,beta)
% Given input of Ho (deepwater wave height), Lo (deepwater wavelength), and
% beta (beach slope, defined as over +/-2sigma about the mean runup line),
% returns the Stockdon 2006 parameterization for runup.
% If the beach is dissipative (iribarren <0.3), will return a special
% dissipative parameterization.
% 
% See "Empirical parameterization of setup, swash, and runup." Stockdon H,
% Holman R, Howd P, Sallenger A, Coastal Engineering 2006 vol: 53 (7) pp:
% 573-588


sqrtHoLo = sqrt(Ho.*Lo);
R2 = 1.1*(0.35*beta.*sqrtHoLo + 0.5*sqrt((0.563*beta.^2 +0.004).*(Ho.*Lo)));% Eq. 19
Sig = 0.06.*sqrtHoLo; % Eq 12
Sinc = 0.75*beta.*sqrtHoLo; % Eq 11
eta = 0.35*beta.*sqrtHoLo; % Eq 10
S = nan(size(Ho)); % preallocate this for dissipative beaches only

%% for Dissipative Beaches
xi = beta./(Ho./Lo).^0.5;
d = find(xi<0.3);

R2(d) = 0.043*sqrtHoLo(d); % Eq 18
Sig(d) = nan;
Sinc(d) = nan;
S(d) = 0.046*sqrtHoLo(d); % Eq 17
