function [k, omega] = dispersion(h,freq)
%% Syntax
% 
% Solves the dispersion relation using the Eckart approximation as a first 
% pass than the Newton-Raphson method up to the 5th decimal place - 
% dispersioncheck: assess accuracy of dispersion solver by plugging 
% wavenumber back into disp relation -> good to 10^-16
%
% [k, omega] = dispersion(h,freq)
% 
% 
% INPUT
%     * h:
%         - type: float
%         - size: h
%         - dimension: meters
%       * freq:
%         - type: float
%         - size: [1 x Nf]
%         - dimension: Hz
%         
%     
%  OUTPUT
%  k: wavenumber
%         - units: 1/m
%  omega: radial frequency values
%         - units: radians
%     
%% 
%% Author
% Athina Lange, SIO July 2019

%% InputParameters
if sign(h) == -1
       h = -h;
end

omega = 2*pi.*freq; % radial frequency
g = 9.81; %m/s

n = 1;
gamma = omega.^2 .*h/g;

%% Doing a first approximation for kh from the Eckart method (1951) 'Surface Waves on Water of Variable Depth'
% omega^2 = gh*sqrt(tanh(omega^2*h/g))
% can also just try placing all beginning values at 0.001

x(:,n) = gamma ./(sqrt(tanh(gamma)));

%% Using a Newton-Raphson approximation
% x_n+1 = x_n - f(x_n)/f'(x_n) where f(x_n) is a tangent line to the
% function at the point x_n in this case 
%epsilon(kh) = omega^2_h/g - (kh)tanh(kh)
% careful with this method when point of interest near local max/min or at
% pts of inflation (good to check with primilarly Eckart result)

x(:,n+1) = x(:,n) - (gamma'-x(:,n).*tanh(x(:,n)))./(-x(:,n).*sech(x(:,n)).^2-tanh(x(:,n)));
for i = 1:size(x,1)
    n = 1;
    while abs(x(i,n+1)-x(i,n)) > 0.000001
        x(i,n+2) = x(i,n+1) - (gamma(i)-x(i,n+1).*tanh(x(i,n+1)))./(-x(i,n+1).*sech(x(i,n+1)).^2-tanh(x(i,n+1)));
        n = n+1;
    end
end

for i = 1:size(x,1)
    aa = x(i,:);
    a(i) = aa(find(aa,1,'last'));
end
k(:) = a./h;
end
