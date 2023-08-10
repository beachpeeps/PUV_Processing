clear all
imgfoldername = ' ̃/Documents/LIDARPressure/Figures';
%% Make some data
Fs = 2;
M = 64;
N = 128;
L = M*N;
t = (0:N-1)*1/Fs;
fN = Fs/2; % Nyquist frequency

%1. Form M sets of data records of length N

omega_b = 0.22*fN;
omega_c = 0.375*fN;
omega_d = omega_b+omega_c;

theta_b = 3*randn(1);
theta_c = 3*randn(1);

testcase = 1; % change this for different tests
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
        case 4
            xi(:,i) = 5*sin(2*pi*omega_c*t);
    end
end
%% 
xi = eta_level(1580).etaAa;
N = 7200/8
M = 7200/N
Fs = 2
fN = Fs/2; % Nyquist frequency

xi = reshape(xi, N, M);

%%

% 2. subtract mean from each record
xi = xi-ones(N,1)*mean(xi);

% 3. apply window to each record
w = hann(N);
w = w./sqrt(mean(w.^2)); % rescale window so that mean(w.ˆ2) is 1

% for m=1:M
% xw(:,m) = xi(:,m).*w;
% end

xw = xi;

% 4. compute the fourier amplitudes
X = fft(xw)/N;
X = X(2:end,:); %ignore the DC component (why?? but okay...)

dt = 0.5; % in cps
df = 1/N;
fny = 1/(2*1);
frequency = [0:df:fny -(fny-df):df:-df];
X(frequency > 0.25) = NaN;
X(frequency < -0.25) = NaN;
%%
% 5. estimate the bispectrum
Bhat = nan(N/4,N/2);
de1 = nan(N/4,N/2);
de2 = nan(N/4,N/2);
bhat = nan(N/4,N/2);

for l=1:N/2
    for k=1:N/4
        if k+l <= N/2 && k <= l
            Bhat(k,l) = 1/M*sum(X(k,:).*X(l,:).*X(k+l,:),2);
            de1(k,l) = sqrt(1/M*sum(abs(X(k,:).*X(l,:)).^2,2));
            de2(k,l) = sqrt(1/M*sum(abs(X(k+l,:)).^2,2));
            bhat(k,l) = abs(Bhat(k,l))/(de1(k,l)*de2(k,l));
        end
    end
end
b2hat = bhat.^2;

% Define frequency axis
df = Fs/N;
% f = Fs/2*linspace(0,1,N/2+1);
f = df:df:Fs/2;

% sum over all segments
Sxxhat = sum(abs(X.^2),2)/M;
Sxxhat = 2*Sxxhat(1:N/2)/df;

figure
semilogy(f,Sxxhat,'k-o')
hold on
plot(PUV_process(1580).Spec.fm, PUV_process(1580).Spec.SSE)
xlabel('Frequency (f)')
ylabel('PSD')

% set figure size
set(gcf, 'units', 'inches', 'pos', [0 0 5 4])
set(gcf, 'PaperPositionMode','auto')
set(gca, 'LooseInset', [0,0,0,0]);

% print to file
% filename = sprintf('PSD case%.0f longN',testcase);
% saveas(gcf,[imgfoldername '/bispectra/' filename '.eps'])

figure
mesh(f,f(1:k),b2hat)
view(22,31)
zlabel('bˆ2(k,l)')
xlabel('Frequency (f)')
zlim([0 1])

% set figure size
set(gcf, 'units', 'inches', 'pos', [0 0 5 4])
set(gcf, 'PaperPositionMode','auto')
set(gca, 'LooseInset', [0,0,0,0]);

% print to file
% filename = sprintf('bispectra case%.0f longN',testcase);
