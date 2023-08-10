clear all
imgfoldername = ' ̃/Documents/LIDARPressure/Figures';

%% Make some data
Fs = 2;
N = 64*128;
t = (0:N-1)*1/Fs;
fN = Fs/2; % Nyquist frequency

omega_b = 0.22*fN;
omega_c = 0.375*fN;
omega_d = omega_b+omega_c;

theta_b = abs(randn(1));
theta_c = abs(randn(1));

testcase = 3; % change this for different tests

switch testcase
    case 1
        theta_d = theta_b+theta_c; % phase−linked theta
        x = cos(2*pi*omega_b*t + theta_b) + ...
            cos(2*pi*omega_c*t + theta_c) + ...
            cos(2*pi*omega_d*t + theta_d) + 0.22*randn(size(t));
    case 2
        theta_d = randn(1); % random theta
        x = cos(2*pi*omega_b*t + theta_b) + ...
            cos(2*pi*omega_c*t + theta_c) + ...
            cos(2*pi*omega_d*t + theta_d) + 0.22*randn(size(t));
    case 3
        theta_d = randn(1); % random theta
        x = cos(2*pi*omega_b*t + theta_b) + ...
            cos(2*pi*omega_c*t + theta_c) + ...
            0.1*cos(2*pi*omega_d*t + theta_d) + ...
            cos(2*pi*omega_b*t + theta_b).*cos(2*pi*omega_c*t + theta_c)...
            + 0.22*randn(size(t));
end



%%
clear all
imgfoldername = ' ̃/Documents/LIDARPressure/Figures';
%% Make some data
Fs = 2;
M = 250; %number of ensembles
N = 128; %length of frame
L = M*N; %total length
t = (0:N-1)*1/Fs;
fN = Fs/2; % Nyquist frequency

% 1. Form M sets of data records of length N

omega_b = 0.22*fN;
omega_c = 0.375*fN;
omega_d = omega_b+omega_c;

theta_b = 3*randn(1);
theta_c = 3*randn(1);

testcase = 3; % change this for different tests
for i=1:M
    switch testcase
        case 1
            theta_d = theta_b+theta_c; % phase−linked theta
            xi(:,i) = cos(2*pi*omega_b*t + theta_b) + ...
                cos(2*pi*omega_c*t + theta_c) + ...
                cos(2*pi*omega_d*t + theta_d) + 0.22*randn(size(t));
        case 2
            theta_d = 3*randn(1); % random theta
            xi(:,i) = cos(2*pi*omega_b*t + theta_b) + ...
                cos(2*pi*omega_c*t + theta_c) + ...
                cos(2*pi*omega_d*t + theta_d) + 0.22*randn(size(t));
        case 3
            theta_d = 3*randn(1); % random theta
            xi(:,i) = cos(2*pi*omega_b*t + theta_b) + ...
                cos(2*pi*omega_c*t + theta_c) + ...
                0.5*cos(2*pi*omega_d*t + theta_d) + ...
                cos(2*pi*omega_b*t + theta_b).*cos(2*pi*omega_c*t + theta_c)...
                + 0.22*randn(size(t));
    end
end


%%

% subtract mean from each record
xi = xi-ones(N,1)*mean(xi);

% apply window to each record
w = hann(N);

for m=1:M
    xw(:,m) = xi(:,m).*w;
end

% normalize the data after windowing it to the variance of original
% timeseries
xw = xw/sqrt(var(xw(:))/var(xi(:)));

% compute the fourier amplitudes
X = fft(xw);
X = X(2:end,:); %ignore the DC component (why?? but okay...)

% Define frequency axis
df = Fs/N;
% f = Fs/2*linspace(0,1,N/2+1);
f = df:df:Fs/2;

% sum over all segments
% Energy Density, units: mˆ2/Hz (Avg Energy per frequency band)
Sxxhat = sum(abs(X.^2),2)/L;
% estimate the bispectrum
Bhat = nan(N/4,N/2);
de1 = nan(N/4,N/2);
de2 = nan(N/4,N/2);
bhat = nan(N/4,N/2);

for l=1:N/2
    for k=1:N/4
        if k+l <= N/2 && k <= l
            Bhat(k,l) = 1/L*sum(X(k,:).*X(l,:).*X(k+l,:),2); %skewness per frequency
            de1(k,l) = Sxxhat(k)*Sxxhat(l)*Sxxhat(k+l); %variance of skewness
            bhat(k,l) = Bhat(k,l)/sqrt(de1(k,l)); %normalized bispectrum
        end
    end
end

% bicoherence squared
b2hat = 2*abs(bhat).^2/M; % Only looking at positive values? I guess.


Sxxhat = df*Sxxhat(1:N/2); % variance preserving
figure
semilogy(f/fN,Sxxhat,'−ok')
xlabel('Normalized Frequency (f/f N)')
ylabel('Energy Density (mˆ2/Hz)')

% set figure size
set(gcf, 'units', 'inches', 'pos', [0 0 5 4])
set(gcf, 'PaperPositionMode','auto')
set(gca, 'LooseInset', [0,0,0,0]);

% % print to file
% filename = sprintf('PSD case%.0f longN',testcase);
% saveas(gcf,[imgfoldername '/bispectra/' filename '.eps'])

figure
mesh(f/fN,f(1:k)/fN,b2hat)
view(22,31)
zlabel('| \Gamma|ˆ2(k,l)')
xlabel('Frequency (f/f N)')
% zlim([0 1])

% set figure size
set(gcf, 'units', 'inches', 'pos', [0 0 5 4])
set(gcf, 'PaperPositionMode','auto')
set(gca, 'LooseInset', [0,0,0,0]);
%
% print to file
filename = sprintf('bispectra case%.0f longN HW',testcase);
saveas(gcf,[imgfoldername '/bispectra/' filename '.eps'])