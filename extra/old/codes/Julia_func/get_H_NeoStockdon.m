function [hsig,hig,eta,higlo,highi,x, depth,CI_hsig,CI_hig] = get_H_NeoStockdon(folderdir,nfft)


% res = load(folderdir,'sur','grd');
% grd = load([folderdir 'grd.mat']);
% sur = load([folderdir 'sur.mat']);

[sur,~] = swash_loadMatFile([folderdir 'sur.mat']);
[grd,~] = swash_loadMatFile([folderdir 'grd.mat']);


bot = -grd.Botlev(1,:);
wlev = squeeze(sur.Watlev(1,:,:));
wlev = wlev-bot'*ones(1,length(wlev));
x = grd.Xp(1,:);
dt = sur.time(2)-sur.time(1);
depth = mean(wlev,2);

trunc50 = 50*60/dt;
trunc50 = round(trunc50);

minutestime = length(wlev(:,end-trunc50:end))./60*dt;
disp(['Length of record is ' num2str(minutestime,'%2.1f') ' minutes'])

[hsig,hig,eta,CI_hsig, CI_hig,higlo,highi] = get_hsig_from_watlev(wlev(:,end-trunc50:end),dt,nfft);


nans = find(isnan(hig));
x(nans) = [];
hsig(nans) = [];
hig(nans) = [];
eta(nans) = [];
higlo(nans) = [];
highi(nans) = [];
CI_hig(nans,:) = [];
CI_hsig(nans,:) = [];


end
