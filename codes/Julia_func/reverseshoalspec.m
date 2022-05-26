function [Ho, spec1Dshoal] = reverseshoalspec(spec1D,h,f,bandwidth,shoaldepth)
% Ho = reverseshoalspec(spec1D,h,f)
% this function reverse shoals the swell frequencies of a spectrum to deep
% water, or shoaling depth of choice.
%
% INPUT:
% spec1D: A 1D spectrum E(time obs,frequency) 
% h: depth of input spectrum
% f: frequency (same size as f used in spec1D)
% shoaldepth: optional shoaling depth
% bandwidth: for use with MOPS data, bandwidth is NOT just f(2)-f(1)!!
%
% OUTPUT:
% Ho: deepwater wave height Ho(time obs)
% spec1Dshoal: reverse-shoaled spectrum
%
% 2019 Julia Fiedler, jfiedler@ucsd.edu

disp('hey')

g = 9.8;

df = bandwidth;
k = get_wavenumber(2*pi*f(:),h);
Cg = get_cg(k(:),h);

revshoal = ones(size(df));

ind = find(f<=0.25); % only want to reverse shoal the swell frequencies

if nargin<5
    revshoal(ind) = 2*Cg(ind)./sqrt(g./k(ind));
elseif nargin == 5
    CgD = get_cg(k(:),shoaldepth);
    revshoal(ind) = Cg(ind)./CgD(ind);
end

spec1Dshoal = bsxfun(@times,revshoal',spec1D);
% Ho = 4*sqrt(nansum(spec1Dshoal(:,ind)* df(ind),2 ) )';
Ho = 4*sqrt(nansum(spec1Dshoal* df,2 ) )';



function cg = get_cg(k,h)

g = 9.81;

% check to make sure k & h are the same size

[nk,mk] = size(k);
[nh,mh] = size(h);

if (nh==1)
   if (mh==1)
      cg = 0.5*(g*tanh(k*h)+g*(k.*h).*(sech(k.*h).^2))./sqrt(g*k.*tanh(k.*h));
      return;
   end;
end;


if (nk ~= nh)
   disp('** Error k & h have wrong sizes');
   return;
end;

if (mk ~= mh)
   disp('** Error k & h have wrong sizes');
   return;
end;

cg = 0.5*(g*tanh(k.*h)+g*(k.*h).*(sech(k.*h).^2))./sqrt(g*k.*tanh(k.*h));

function k = get_wavenumber(omega,h)

% returns the wavenumber of the gravity wave
% dispersion relation, by using newtons method

% the initial guess will be the shallow water wavenumber

g = 9.81;

k = omega./sqrt(g*h);


f = g*k.*tanh(k.*h) - omega.^2;

while max(abs(f))>1e-10
  dfdk = g*k.*h.*(sech(k.*h)).^2 + g*tanh(k.*h);
  k = k - f./dfdk;
  f = g*k.*tanh(k.*h) - omega.^2;
end



