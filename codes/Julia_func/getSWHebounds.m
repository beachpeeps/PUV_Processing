function [ Hs, lb, ub, edof ] = getSWHebounds( e, dof, q, df )
% [ Hs, lb, ub, edof ] = getSWHebounds( e, dof )
% This function computes the effective degrees of freedom (edof) from a
% spectrum, following Elgar 1987.
%
% Inputs:
%   e - Energy(time,freq), units: m^2
%   dof - degree of freedom per freq-band
%   q - percentile to determine ub,lb (e.g. 0.90 or 0.68)
%   df - frequency bin size
%
% Outputs:
%   lb - lower Hs bound
%   ub - upper Hs bound
%   Hs - Hsig
%   edof - Effective dof used
%
% References:
%   ELGAR 1987, for details of un-biased EDOF estimate
%   used in Fiedler et al. 2018 Coastal Engineering, section 4.5
%

% copyright 2017 Julia W Fiedler jfiedler@ucsd.edu

% percent on either side
a = (1-q)/2;

% Get Hsig
Hs = 4*sqrt(sum(e.*df));

% Determine edof                             
edof = dof*sum(e).^2./sum(e.^2);    %estimate effective DOF
edof = edof/(1+2/dof);                    %UNBIAS THE ESTIMATE 

% Generate normalized ub,lb

lb = edof/chi2inv(1-a,edof);
ub = edof/chi2inv(a,edof);

% Multiply by Hsig
lb = sqrt(lb).*Hs;
ub = sqrt(ub).*Hs;

% % Multiply by Hsig
% lb = lb.*Hs;
% ub = ub.*Hs;

end