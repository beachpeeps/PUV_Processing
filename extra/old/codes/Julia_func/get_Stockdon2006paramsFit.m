function [R2mdl,Sigmdl,Sincmdl,Etamdl,sqrtHoLo,mdlSig,mdlSinc,mdlEta] = get_Stockdon2006paramsFit(Ho,Lo,beta,Sig,Sinc,etaRunup)
% [R2,Sig,Sinc,eta,sqrtHoLo] = get_Stockdon2006params(Ho,Lo,beta)
% Given input of Ho (deepwater wave height), Lo (deepwater wavelength), and
% beta (beach slope, defined as over +/-2sigma about the mean runup line),
% returns the tuned Stockdon 2006 parameterization for runup.
% 
% See "Empirical parameterization of setup, swash, and runup." Stockdon H,
% Holman R, Howd P, Sallenger A, Coastal Engineering 2006 vol: 53 (7) pp:
% 573-588


sqrtHoLo = sqrt(Ho.*Lo);

mdlSinc = fitlm(beta.*sqrtHoLo./2,Sinc./2,'Intercept',false);
mmSinc = mdlSinc.Coefficients.Estimate;
Sincmdl = mmSinc*beta.*sqrtHoLo ;

mdlEta = fitlm(beta.*sqrtHoLo,etaRunup,'Intercept',false);
mmEta = mdlEta.Coefficients.Estimate;
Etamdl = mmEta*beta.*sqrtHoLo;


mdlSig = fitlm(sqrtHoLo./2,Sig./2,'Intercept',false);
mmSig = mdlSig.Coefficients.Estimate;
Sigmdl = mmSig*sqrtHoLo;


R2mdl = 1.1*(Etamdl + 0.5*sqrt(Sigmdl.^2+Sincmdl.^2));
% R2 = 1.1*(0.35*beta.*sqrtHoLo + 0.5*sqrt((0.563*beta.^2 +0.004).*(Ho.*Lo)));% Eq. 19
% Sig = 0.06.*sqrtHoLo; % Eq 12
% Sinc = 0.75*beta.*sqrtHoLo; % Eq 11
% eta = 0.35*beta.*sqrtHoLo; % Eq 10
% S = nan(size(Ho)); % preallocate this for dissipative beaches only


