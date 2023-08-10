ffunction [etaAa, etaAss, etaAig, frequency, time]= etasolver_detide(pm, tm, interval, h, varargin)
%% Syntax
% Calculating the sea surface elevation from bottom-mounted pressure
% sensors (assuming no current) for a SS band [1/25 - 0.35Hz] and an IG band
% [1/interval - 1/25Hz]
% removes S1 and M2 tides as well as mean sea level in that 1.5h interval
% computes wavenumber from dispersion.m (Eckart and Newton-Raphson method)
%
% Requires: dispersion.m
% 
%
% [etaAa, etaAss, etaAig, frequency, time]= etasolver(pm, tm, 5400)
% [etaAa, etaAss, etaAig, frequency, time]= etasolver(pm, tm, 5400, 1619)
% 
%
% INPUT
%       * pm: pressure sensor data (bottom-mounted sensor)
%               - type: float
%               - size: [Nlength x Nseg] : 
%                           Nlength (size of interval), Nseg (number of segments)
%               - dimension: m [p/(rho*g)]
%       * tm: time vector
%               - type: datenum
%               - size: [Nlength x Nseg]
%               - dimension: s 
%       * interval: size of interval input (5400 or 1800s)
%       * varargin: 
%               - sid: which 1.5h segment to be plotted -> default 1
% 
% OUTPUT
%  etaAa: sea surface elevation for SS and IG band (1/5400 - 0.35Hz)
%           - units: m
%  etaAss: sea surface elevation for SS (1/25 - 0.35Hz)
%           - units: m
%  etaAig: sea surface elevation for IG band (1/5400 - 1/25Hz)
%           - units: m
%  frequency : frequency vector [0 :df: f_ny  -(f_ny-df) : df:-df]
%       - units: Hz
%  time : time vector
%       - units: s
%
%% 
%% Author
% Athina Lange, SIO July 2019
%% 

p = inputParser();
p.CaseSensitive = false;
p.addOptional('sid',1);
p.parse(varargin{:});

sid = p.Results.sid;

time = tm;%datetime(tm, 'ConvertFrom', 'datenum');
%interval = 5400; 
[m,n] = size(pm);
timedat = (1:m)';
%% Remove mean water level and S1 and M2 tides from pressure data by least squares fitting
% See DATA1 notes
%h = mean(pm, 1);
% pd = pm% - h; already detided;
% x = NaN(5,n);
% trend = NaN(m,n);
% A = [ones(length(timedat),1), ...
%             sin(2*pi*timedat/24/3600),cos(2*pi*timedat/24/3600),... %S1
%             sin(2*pi*timedat/12.42/3600),cos(2*pi*timedat/12.42/3600)]; %M2
% for i = 1:size(pm,2)
%     x(:,i) = inv(A.' * A) * A.' * pd(:,i);
%     trend(:,i) = A * x(:,i);
% end

pd = pm;%pd - trend;
pd=reshape(detrend(pd(1:floor(m/interval)*interval,:)),interval,n*floor(m/interval));
pms = reshape(pm(1:floor(m/interval)*interval,:), interval, n*floor(m/interval));
times=reshape(time(1:floor(m/interval)*interval,:),interval,n*floor(m/interval));

%h = mean(pms,1); % mean water level = 7.7286m ~ 8m
%% Frequency and wavenumber
dt = nanmean(median(diff(time)),2); % in cps
% df = 1/size(pms,1);
% fny = 1/(2*1);
% frequency = [0:df:fny -(fny-df):df:-df];

nfft=size(pms,1);
df = 2/(nfft-1);   % This is the frequency resolution
nnyq = nfft/2 +1;
frequency = [0:df:(nnyq-1)*df -(nnyq-2)*df:df:-df];


k = NaN(size(pd,1), size(pd,2));
for i = 1:size(h,2)
    [k(:,i), ~] = dispersion(h(i),frequency);
    %[~, k(:,i)] = dispr(h(i), frequency);
end

%% Pressure response factor (inverse)
flim = [0.004 0.04 0.25]; % limits of IG and SS bands
igl = find(abs(frequency) <= flim(1));
igh = abs(frequency) > flim(2);
ssl = abs(frequency) <= flim(2);
ssh = find(abs(frequency) >= flim(3)); %same

% 
Kpi = cosh(k(:,:).*h(:)');
% All frequencies
Kpi(igl,:) = 0;
Kpi(ssh,:) = 0;
% SS
KpiAss = Kpi;
KpiAss(ssl,:) = 0;
KpiAss(ssh,:) = 0;
% IG
KpiAig = Kpi;
KpiAig(igh,:) = 0;
KpiAig(igl,:) = 0;

%% Calculating eta
pdF = fft(pd);

% All
etaF = pdF.*Kpi;
etaA = real(ifft(etaF));
etaAa = etaA;%reshape(etaA, interval*floor(m/interval), n);

%SS
etaFss = pdF(:,:).*KpiAss(:,:);
etaAss = real(ifft(etaFss));
%etaAss = reshape(etaAss, interval*floor(m/interval), n);

%IG
etaFig = pdF(:,:).*KpiAig(:,:);
etaAig = real(ifft(etaFig));
%etaAig = reshape(etaAig, interval*floor(m/interval), n);

time = reshape(times, interval*floor(m/interval),n);

end