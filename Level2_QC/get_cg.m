% Falk Feddersen (c) 2001
%
% function that takes the vector wavenumbers and depths
% and returns the group velocity
%
% function cg = funwaveC_get_cg(k,h)

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