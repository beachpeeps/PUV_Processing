function [siglevel] = significance(xlength,Fs,fmax,ensemblelength)
%finds the 95% significance level by brute force

ind=1;
siglevel = nan(100,1);
while ind<=100

% generate some white noise, unit variance
x = randn(xlength,1);

% bicoherence of data
[~,~,~,~,b2hat,~,~,~,~,~,~] = bispectra_juliatrial(x,Fs,fmax,ensemblelength);

fivepercent = round(0.5*sum(~isnan(b2hat(:))));
b2hatsort = sort(b2hat(~isnan(b2hat)),'descend');

siglevel(ind) = b2hatsort(fivepercent);
ind = ind+1;
end
siglevel = mean(siglevel);