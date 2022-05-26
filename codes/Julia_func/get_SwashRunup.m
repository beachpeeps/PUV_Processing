function get_SwashRunup(filedir,casenum,threshold)
% folder = '~/Documents/SWASH_runup/agate_21_10layer_loglaw0008_numerics2/';

if nargin<1
    filedir = 'agate_22/';
    casenum = 22;
    threshold = 0.1;
end

folderdir = '~/Documents/SWASH_runup/models/';
% folder = [folderdir filedir];
%     '~/Documents/SWASH_runup/agate_22/';
% filedir = 'agate_21_10layer_loglaw0008_numerics2/';

filename = 'agate.mat';

% casenum = 22;

[Spec,Info,Bulk,Tseries] = get_runup_agate(folderdir,filedir,filename,casenum,threshold);

Rmin = min(Tseries.Xrunup);
Rxmean = mean(Tseries.Xrunup);

thr = sprintf('_%02.0fcm',threshold*100);
modnum = filedir(end-3:end-1);
savename = ['SwashRunup_' num2str(casenum)];

save(['~/Documents/SWASH_runup/mat/' savename modnum thr ])

disp(['Swash runup stats saved to' savename modnum thr])
figure
plot(Tseries.Xrunup,Tseries.T)