function[f]=fourier(varargin)
%FOURIER  The one-sided Fourier frequencies for a given length time series.
%
%   F=FOURIER(N) returns the one-sided (or positive) Fourier frequencies 
%   for a time series of length N.  
%
%   F is a radian or angular frequency so that the Nyquist is at pi.
%
%   F=FOURIER(DT,N) uses sample rate DT in calculating the frequencies, so
%   that the Nyquist will be at pi/DT.
% 
%   Note that the highest resolvable frequency, MAX(F), differs for even or
%   odd N.  For even N, it is the Nyquist pi/DT, but for odd N the Nyquist
%   is not resolved and the frequency is (N-1)/N * pi/DT.
%
%   F is of length N/2+1 for even N, and (N+1)/2 for odd N. 
%
%   Usage: f=fourier(N);
%          f=fourier(dt,N);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2011--2012 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(varargin{1}, '--t')
    fourier_test,return
end

if nargin==1
    N=varargin{1};
    dt=1;
else
    dt=varargin{1};
    N=varargin{2};
end

if iseven(N)
    f=(0:1./N:1/2)';
elseif isodd(N)
    f=(0:1./N:1/2*frac(N-1,N))';
end

f=pi*f./dt;

function[]=fourier_test

reporttest('FOURIER is length N/2+1 for even N',length(fourier(10))==6);
reporttest('FOURIER is length (N+1)/2 for odd N',length(fourier(11))==6);

%Even and odd length frequencies behave differently, which you can see with 
%abs(fft(cos(pi*[1:10]'))).^2
%abs(fft(cos(pi*[1:11]'))).^2
%abs(fft(cos((10/11)*pi*[1:11]'))).^2
