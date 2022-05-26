function [hsig,hig,eta,higlo,highi,x, datahour, bulkinfo,depth,CI_hsig,CI_hig] = get_H_spatial(folderdir,casenum,method,type,nfft,test)




switch type
    case 'ebalance'
        switch test
            case 'variable'
                res = load([folderdir 'agate.mat'],'sur','grd');
            case 'dfdx'
                res = load(folderdir,'sur','grd');
        end
        bot = -res.grd.Botlev(1,:);
        wlev = squeeze(res.sur.Watlev(1,:,:));
        wlev = wlev-bot'*ones(1,length(wlev));
        x = res.grd.Xp(1,:);
        dt = res.sur.time(2)-res.sur.time(1);
        depth = mean(wlev,2);
    case 'regular'
        [res,~] = swash_loadMatFile([folderdir 'agate.mat']);
        bot = -res.Botlev(1,:);
        wlev = squeeze(res.Watlev(1,:,:));
        wlev = wlev-bot'*ones(1,length(wlev));
        x = res.Xp(1,:);
        dt = res.time(2)-res.time(1);
        depth = mean(wlev,2);
end





trunc50 = 50*60/dt;
trunc50 = round(trunc50);

minutestime = length(wlev(:,end-trunc50:end))./60*dt;
disp(['Length of record is ' num2str(minutestime,'%2.1f') ' minutes'])

switch method
    case {'shoreward'}
        U = squeeze(res.vel_x(1,:,:));
        [hsig,hig,eta,~,~,higlo,highi] = get_hsigshoreward_from_watlev(wlev(:,end-trunc50:end),U,depth,dt);
    case {'total'}
        %             [hsig,hig,eta,CI_hsig, CI_hig,higlo,highi]
        [hsig,hig,eta,CI_hsig, CI_hig,higlo,highi] = get_hsig_from_watlev(wlev(:,end-trunc50:end),dt,nfft);
end

nans = find(isnan(hig));
x(nans) = [];
hsig(nans) = [];
hig(nans) = [];
eta(nans) = [];
higlo(nans) = [];
highi(nans) = [];
CI_hig(nans,:) = [];
CI_hsig(nans,:) = [];


switch type
    case 'ebalance'
        switch test
            case 'variable'
                filename = sprintf('agate%02.0f.sws',casenum);
                fid = fopen([folderdir filename]);
                for i=1:5
                    fgetl(fid);
                end
                datahour = datenum(fgetl(fid),'$ DD-mmm-yyyy HH:MM');
                bulkinfo = fgetl(fid); bulkinfo = bulkinfo(3:end);
            case 'dfdx'
                datahour = nan;
                bulkinfo = nan;
        end
    case 'regular'
        filename = sprintf('agate%02.0f.sws',casenum);
        fid = fopen([folderdir filename]);
        for i=1:5
            fgetl(fid);
        end
        datahour = datenum(fgetl(fid),'$ DD-mmm-yyyy HH:MM');
        bulkinfo = fgetl(fid); bulkinfo = bulkinfo(3:end);
end
end
