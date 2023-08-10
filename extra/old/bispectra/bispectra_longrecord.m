%% bispectra longrecord.m
% function bispectra longrecord(nsensor,hr,min,start,Nlength)
clf

sensor = {'P01','P02','P03','P04','P05','P06','P07','P08','P09','P10'};
dt = 0.5; % in seconds
Fs = 1/dt; % sampling rate

n = 10; % sensor select
hr = 3; % number of hours to link together

dir str = 'cd  ̃/Documents/LIDARPressure/data/Processed Data/';
imgfoldername = ' ̃/Documents/LIDARPressure/Figures';

ss = strcat(dir str,sensor{n}); % open folder to get sensor data
eval(ss)

file list = dir('*00.dat'); % lists all 00.dat files in the directory
endtime = 59*60*2; % for truncating the time series
minremain = endtime/2/60; % number of minutes of data
%%

% start = 131; %this index corresponds to the max sig. wave height at
% P10 over the experimental period
Maxb2 = nan(1,300);
start = 1;
% while start<300
start = 32;
presAll = [];

for x=start:start+hr
    [time,pres] = getpres(x, file list);
    [H,presFilt] = pcorrect(n,x,pres,Fs); % corrects for burial
    % Error if pressure sensor is above water level (dry beach)
    % if H<0
    % error(strcat('H should not be less than 0, Time : ',...
    % datestr(start time)))
    % end
    %% Filter data
    % cutoff = 0.31; % beyond this is noise
    % order = 7;
    % if H>=0
    % [time,presFilt] = myButterfilt(time,presFilt,cutoff,order,Fs);
    % end
    %% Use data to make white noise


    presFilt = presFilt(1:endtime); % truncate time series


    % Demean time series prior to padding with zeros, otherwise will
    % get ringing due to box car shape
    presFilt = presFilt-mean(presFilt);
    presF temp = [presFilt; zeros(length(time)-length(presFilt),1)];
    presAll = [presAll; presF temp];


end

% bispectrum of data
fmax = 0.4; % only look at data to 0.4 Hz
Nlength = 512; % ensemble length
[f,fN,fstop,Sxxhat,b2hat,bphase,N,M,varx] = mybispectra(presAll,Fs,fmax,Nlength);
siglevel = significance(length(presAll),Fs,fmax,Nlength);


%%
% figure
subplot(3,1,1)
semilogy(f(1:fstop),Sxxhat(1:fstop),'k')
ylim([min(Sxxhat(1:fstop)) max(Sxxhat(1:fstop))])

title(sprintf('Sensor %.0f',n))
xlabel('Frequency (Hz)')
ylabel('Energy (mˆ2/Hz)')

dof = 2*M; % 2 times the number of ensembles
alpha = 0.05;
chi2 = [chi2inv(alpha/2,dof) chi2inv(1-alpha/2,dof)];
CI_chi2 = [(1/chi2(1)).*dof 1/chi2(2).*dof];
hold on
CI_chi2 = CI_chi2*median(Sxxhat(1:fstop));

semilogy([f(5) f(5)],CI_chi2,'k','LineWidth',2)



subplot(3,1,2:3)
b2hat(b2hat< siglevel) = NaN;

pcolor(f(1:fstop),f(1:fstop/2),b2hat)
shading flat
xlabel('Frequency (Hz)')
ylabel('Frequency (Hz)')
caxis([0 1])
h = colorbar('Location','SouthOutside');
xlabel(h,'bˆ2(k,l)')

% set figure size
% set(gcf, 'units', 'inches', 'pos', [0 0 5 10])
% set(gcf, 'PaperPositionMode','auto')
% set(gca, 'LooseInset', [0,0,0,0]);
% set(gcf, 'Renderer', 'opengl')

% % print to file
filename = sprintf('bicoh P%02.0f',n);
saveas(gcf,[imgfoldername '/bispectra/' filename '.eps'],'psc2')