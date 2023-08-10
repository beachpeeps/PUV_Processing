function [mem] = get_mem(a1,b1,a2,b2)
% GETMEM returns the varibale array "mem", a N by 360 matrix of normalized
% MEM directional spectra for the input a1(N)...b2(N) normalized 
% directional fourier coefficients where N is the number of frequency 
% bands.  To make a complete 2d freq-directional spectrum with energy
% density of m^2/deg-Hz you must multiply the N normalized distributions
% by the wave energy at each frequency.
%
%  [mem] = getmem(a1,b1,a2,b2)
%  [mem] = getmem(d)  
%
% Maximum entropy method (MEM) for the estimation of a directional
% distribution from 1st and 2nd moment normalized directional fourier
% coefficients.
%
% Lygre and Krogstad, JPO v16 1986: NOTE - there is a typo in the
% paper...BOTH exponetials in eq. 13 should have negative signs.
% This is fixed in a later in paper by Krogstad in a wave spectra
% conference proceedings , Tom Herbers has the book)
%

% Check number of inputs.
if nargin > 4
    error('getmem:TooManyInputs', ...
        'requires at most 4 inputs');
end

%... switch to Lygre & Krogstad notation
switch nargin
    case 1
        d1 = a1(1,:);
        d2 = a1(2,:);
        d3 = a1(3,:);
        d4 = a1(4,:);
    case 4
        d1=a1;
        d2=b1;
        d3=a2;
        d4=b2;
end



c1=complex(1.,0).*d1+complex(0,1.).*d2;
c2=complex(1.,0).*d3+complex(0,1.).*d4;

p1=(c1-c2.*conj(c1))./(1-abs(c1).^2);
p2=c2-c1.*p1;

x=1.-p1.*conj(c1)-p2.*conj(c2);

% calculate MEM using 1 deg resolution
a=1:360;
a=a.*pi/180;
e1=complex(1.,0)*cos(a)-complex(0,1.)*sin(a);
e2=complex(1.,0)*cos(2*a)-complex(0,1.)*sin(2*a);
y=abs(complex(1.,0)-p1'*e1-p2'*e2).^2;

mem=abs((x'*ones(1,360))./y);

% normalize MEM distribution
mem=((1./sum(mem'))'*ones(1,360)).*mem;

%  direction indices are reveresed after matrix
%  calcs, so flip directions back to native 
%  direction coordinates.
mem=fliplr(mem);

return
