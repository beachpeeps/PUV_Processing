% need Grid, Sur, Tseries, Info
clear

% load('~/Documents/SWASH_runup/mat/processed/runupStats', 'Tseries','Info','beta')
load('~/Documents/NeoStockdon/processed/bulkrunup.mat', 'Tseries','Info','Bulk')

savedir = '~/GoogleDriveUCSD/MOPS_Climatology/';
savename = 'runupMOVIE';


%% get variables for plotting

for i = 1:2
    grd(i) = load(Info(i).grdDir);
    sur(i) = load(Info(i).surDir);
    
%     xgrid(i).x = read_swash_coord(Info(i).grdDir(1:end-7));
end

%% find locations where dx switches to finer grid
% dx2 = diff(diff(xgrid(5).x));
% dx2ind = find(abs(dx2)>1e-4)+2;
% xdxind = xgrid(5).x(dx2ind);
%% truncate in time for watlev figure
dt = sur(1).time(2)-sur(1).time(1);
dt = floor(dt*10)./10; % want accuracy only to 0.1
%take 1 50 minute chunk for proper spectra comparison
totalSec = floor(50*60/dt);

% truncate all records to 50 minutes
for i=1:2
    sur(i).Watlev = sur(i).Watlev(1,:,end-totalSec-1:end);
end

%%
figwidth = 7;
figheight = 6;
nrow = 2;
ncol = 1;
units = 'inches';
cbar = 0;
[hFig, ax] = makeFig(figwidth,figheight,nrow,ncol,units,cbar);
hFig.Color = 'w';

ax(2).Position(4) = ax(2).Position(4)-0.05;

leftside = ax(1).Position(1)+ax(1).Position(3);
ax(3) = axes('position',[leftside-0.05 ax(1).Position(2) 0.05 ax(1).Position(4)]);
ax(1).Position(3) = ax(1).Position(3) - ax(3).Position(3)-0.01;
ax(3).Box = 'on';
ax(3).YTickLabel = []; ax(3).XTickLabel = [];

cmap = lines(5);


axes(ax(1))
plot(grd(1).Xp,-grd(1).Botlev,'k')
hold on
% plot([xdxind(end) xdxind(end)],[-5 5],':k')
for i=1:2
    ha1(i) = plot(grd(i).Xp,squeeze(sur(i).Watlev(1,:,1)),'color',cmap(i,:));
end

indStart = 1600;
indStop = 3001;
Tvec = Info(1).datahour+seconds(Tseries(1).T)-minutes(10);
axes(ax(3))
hold on
%%
axes(ax(2))
hold on
for i=1:2
    plot(Tvec(indStart:indStop),Tseries(i).Zrunup(indStart:indStop),'color',[cmap(i,:) 0.2])
    hold on
    ha2(i) = plot(Tvec(indStart:indStart),Tseries(i).Zrunup(indStart:indStart),'color',cmap(i,:));
    ha3(i) = plot([0 1],Tseries(i).Zrunup(indStart)*ones(1,2),'color',cmap(i,:),'parent',ax(3));
    ha4(i) = plot([-100 100],Tseries(i).Zrunup(indStart)*ones(1,2),':','color',cmap(i,:),'parent',ax(1));
    
end

ax(1).XLim = [-60 65];
ax(1).YLim = [-2 5];
ax(2).XLim = [Tvec(indStart) Tvec(indStop)];
ax(2).YLim = [-2 4];
ax(3).YLim = ax(1).YLim;
ax(3).XLim = [0 1];
ax(1).YLabel.String = 'Z (m)';
ax(1).XLabel.String = 'X (m)';
ax(2).YLabel.String = 'Z (m)';
ax(1).Title.String = ['Runup, \beta = ' num2str(Bulk(1).beta,'%02.2f')];

% legend
% indOrder = [3 2 1 ]; %because I arranged the runupStats wrong. TODO: FIX THIS

legStr = {'unimodal','bimodal'};
yloc = [0.29 0.22 0.15 0.08 0.01];
yloc = [0.29 0.22];

for i=1:2
    ht = text(0.83,yloc(i),legStr{i},'units','normalized','parent',ax(1),'color',cmap(i,:),'fontsize',14);
end

igif = 1;
for i= indStart:5:indStop
% for i=indStop
    delete(ha1)
    delete(ha2)
    delete(ha3)
    delete(ha4)
    for ii=1:2
%         ha1(ii) = plot(grd(ii).Xp,squeeze(sur(ii).Watlev(1,:,i))-Info(ii,betan).thinLayer,'parent',ax(1),'color',cmap(ii,:));
        ha1(ii) = plot(grd(ii).Xp,squeeze(sur(ii).Watlev(1,:,i)),'parent',ax(1),'color',cmap(ii,:));
        ha2(ii) = plot(Tvec(indStart:i),Tseries(ii).Zrunup(indStart:i),'color',cmap(ii,:),'parent',ax(2));
        ha3(ii) = plot([0 1],Tseries(ii).Zrunup(i)*ones(1,2),'color',cmap(ii,:),'parent',ax(3));
        ha4(ii) = plot([-100 100],Tseries(ii).Zrunup(i)*ones(1,2),':','color',cmap(ii,:),'parent',ax(1));
    end
    
    print(['~/Desktop/temp/Frame ' num2str(igif)], '-dpng','-r150');
    drawnow
    F = getframe(hFig);
    [im{igif},map] = frame2im(F);
    igif = igif+1;
end


% print(hFig, '-dpdf', [savedir savename]);


%%
gifname = [savedir savename '.gif'];
for idx = 1:igif-1
    [A,~] = imread(['~/Desktop/temp/Frame ' num2str(idx) '.png']);
    [X,map] = rgb2ind(A,256);
    if idx == 1
        imwrite(X,map,gifname,'gif','LoopCount',Inf,'DelayTime',0);
    else
        imwrite(X,map,gifname,'gif','WriteMode','append','DelayTime',0);
    end
end



%
%
%
% %%
% wlev10 = Wlev-thinLayer;
% wlev10(wlev10>0.05) = nan;
% clf
% pcolor(grd.Xp,sur.time(5000:6000),wlev10(5000:6000,:))
% shading flat
%
% xlim([-100 0])
% %%
% figure
% subplot(2,1,1)
% h = plot(grd.Xp,wlev10(5000:10000,:),'color',[0 0 0 0.2]);
% xlim([-50 0])
%
% hold on
% plot(grd.Xp,nanmean(wlev10(5000:10000,:)),'r')
% subplot(2,1,2)
%
% percentData = sum(~isnan(Wlev))./length(sur.time);
% ind80percent = find(percentData<0.8,1);
% meanRunupShape = nanmean(Wlev);
% thinLayer = meanRunupShape(ind80percent);
%
% plot(grd.Xp,percentData)
% xlim([-50 0])