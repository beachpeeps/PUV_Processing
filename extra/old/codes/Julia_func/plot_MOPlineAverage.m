clear
load ~/Documents/NeoStockdon/data/D0045_sand_levels.mat
t = datetime(dates,'InputFormat','yyyyMMdd');

% Mean Seal Level in NAVD88
NAVD88mllw = 0.058;
MSLmllw = 0.832;
MSLnavd88 = MSLmllw-NAVD88mllw;
MHWmllw = 1.402;
MHWnavd88 = MHWmllw-NAVD88mllw;

% tNourishment = datetime(2013,2,1):datetime(2012,9,1)+calyears(3);
% 
% indNotNourishment = find(~ismember(t,intersect(t,tNourishment)));
% indNotNourishment(18) = [];
% 
% indNourishment = find(ismember(t,intersect(t,tNourishment)));
% 
smoothwidth = 10;
% 
% tt = t(indNotNourishment);
% D0045zz = D0045z(:,indNotNourishment);
% 
% 
% tN = t(indNourishment);
% D0045zN = D0045z(:,indNourishment);

tt = t;
tt(18) = [];
D0045zz = D0045z;
D0045z(:,18) = [];

monthnum = t.Month;
yearnum = t.Year;
years = unique(yearnum);
cmap = cmocean('phase');

mshift = nan(1,length(t));
% months departure from August, this is beyond klugey but whatever
for i=1:length(t)
    if monthnum(i)>2
        mshift(i)  = abs(monthnum(i)-8);
    elseif monthnum(i)<=2
        mshift(i) = monthnum(i)+4;
    end    
end


clf
% yearind = [3:9];
% ny = numel(yearind);
% ny = numel(unique(mshift));

% cmap = cmocean('thermal',ny);
cmap = cmocean('-balance');

n = 1;
for i= [7 6 5 4 3 2 1]
%     monthind = find(yearnum == years(yearind(n)));
    monthind = find(mshift == i-1);

    plot(D0045crs,D0045zz(:,monthind),'color',[cmap(round(i/7*256),:) 0.1])
    hold on
    n = n+1;
end



% find all profiles that are eroded
%  use the profiles that intersect with MSL before x=0

L1 = [[-350 50];[1 1]*MSLnavd88];
for i=1:length(tt)
sandLevel = [D0045crs; D0045z(:,i)'];
try
temp = InterX(L1,sandLevel);
beachWidth(i) = temp(1);
catch
end

end

eroded = find(beachWidth>10);

% plot(D0045crs,D0045z(:,eroded),'color',[0 0 0 0.2])


% extend to z=-10
erodedProfile = nanmean(D0045z(:,eroded),2);

% take slope from last two points
ind = find(~isnan(erodedProfile),1);
mdl = fitlm(D0045crs(ind:ind+1),erodedProfile(ind:ind+1));
mb = mdl.Coefficients.Estimate;
% 
% newslope = mb(1)+D0045crs*mb(2);
% xind = knnsearch(newslope(:),-10);

%extend to x=-15;
x15 = (-15-mb(1))./mb(2);
x15 = round(x15/100)*100;
%%
newx = x15:5:D0045crs(1)-5;
newxonshore = D0045crs(end)+5:5:150;
D0045crsNew = [newx D0045crs newxonshore];
newslope = mb(1)+D0045crsNew*mb(2);
% erodedProfile = [-15; erodedProfile];
erodedProfile = [nan(size(newx))'; erodedProfile];
ind = find(~isnan(erodedProfile),1);

erodedProfile(1:ind) = newslope(1:ind);

% extend to z=10
ind = find(~isnan(erodedProfile),1,'last');
mdl = fitlm(D0045crsNew(ind-1:ind),erodedProfile(ind-1:ind));
mb = mdl.Coefficients.Estimate;

newslope = mb(1)+D0045crsNew*mb(2);
xind = knnsearch(newslope(:),15+MSLnavd88);

erodedProfile(ind:xind) = newslope(ind:xind);
%%

% smooth the eroded profiles
zSmooth(:,1) = smoothdata(erodedProfile,'gaussian',10,'includenan'); % 50m smoothing
zSmooth(:,2) =  smoothdata(erodedProfile,'gaussian',4,'includenan'); %20m smoothing
zSmooth(:,3) = erodedProfile; 

% chunk out areas for smoothing
% let x=-300 to x = -50 have 50 m smoothing
offshore = find(D0045crsNew>-450 & D0045crsNew<-20);
zSmooth(offshore,2:3) = NaN;
smoothErodedProfile = nanmean(zSmooth,2);

plot(D0045crsNew,smoothErodedProfile,'.k','linewidth',2)


ax = gca;

plot(ax.XLim,[1 1]*MSLnavd88,':k')
plot(ax.XLim,[1 1]*MHWnavd88,':k')

plot([0 0],ax.YLim,':k')
text(0,-6,'MEAN SHORELINE','rotation',90,'background','w')
text(-650,MSLnavd88,'MSL','background','w')
text(-650,MHWnavd88,'MHW','background','w')

xlabel('X (m)')
ylabel('Z (navd88,m)')

text(-650,2.5,...
    ['Representative Eroded Profile ' datestr(tt(1)) ' - '  datestr(tt(end))],...
    'fontsize',16)

figureDir = '~/GoogleDriveUCSD/MOPS_Climatology/';
figureName = [figureDir 'IBprofile3.jpeg'];
% 
% print('-djpeg',figureName,'-r300')

%% save data for SWASH runs
z = smoothErodedProfile;
x = D0045crsNew; 

indStart = find(~isnan(z),1);
indStop = find(~isnan(z),1,'last');

z = z(indStart:indStop);
x = x(indStart:indStop);
xEroded = x;
zEroded = z;

xx = [min(x):2:-202 -200:1:-20 -19:0.5:0 0.25:0.25:max(x)];
zz = interp1(x,z,xx);

bathy.x = xx;
bathy.z = zz;
bathy.dx = 'variable';
bathy.description  = 'MOP line D0045, eroded profile';


% save('~/Documents/NeoStockdon/data/D0045profile_eroded_15.mat','bathy')
%%

% colormap(cmap)
% c = colorbar(gca,'location','northoutside');
% % c.Ticks = linspace(0.5/ny,1-(0.5./ny),ny);
% % c.TickLabels = num2str(years(yearind));
% c.Ticks = [0 1];
% c.TickLabels = {'summer','winter'};

z = nanmedian(D0045zz(:,ismember(monthnum,[ 11 12 1 2 3])),2);
zSmooth = smoothdata(z,'gaussian',smoothwidth);


% winter profile
% plot(D0045crs,zSmooth,'color',cmap(end-1,:),'linewidth',3)
% summer profile
% plot(D0045crs,nanmean(D0045zz(:,intersect(monthnum,[6 7 8 9 10])),2),'color',cmap2(2,:),'linewidth',3)

%% Nourished Profile

nourished = find(beachWidth<0);

% plot(D0045crs,D0045z(:,nourished),'color',[0 0 0 0.2])

% extend to z=-15
nourishedProfile = nanmean(D0045z(:,nourished),2);
plot(D0045crs,nourishedProfile,'r','linewidth',1)


% take slope from last two points
ind = find(~isnan(nourishedProfile),1);
mdl = fitlm(D0045crs(ind:ind+1),nourishedProfile(ind:ind+1));
mb = mdl.Coefficients.Estimate;

newslope = mb(1)+D0045crs*mb(2);
x15 = (-15-mb(1))./mb(2);
% xind = knnsearch(newslope(:),-10);
newx = x15:5:D0045crs(1)-5;
newxonshore = D0045crs(end)+5:5:200;
D0045crsNew = [newx D0045crs newxonshore];
newslope = mb(1)+D0045crsNew*mb(2);
% nourishedProfile(xind:ind) = newslope(xind:ind);

nourishedProfile = [nan(size(newx))'; nourishedProfile];
ind = find(~isnan(nourishedProfile),1);

nourishedProfile(1:ind) = newslope(1:ind);


% extend to z=10
ind = find(~isnan(nourishedProfile),1,'last');
mdl = fitlm(D0045crsNew(ind-1:ind),nourishedProfile(ind-1:ind));
mb = mdl.Coefficients.Estimate;

newslope = mb(1)+D0045crsNew*mb(2);
xind = knnsearch(newslope(:),10+MSLnavd88);

nourishedProfile(ind:xind) = newslope(ind:xind);


% smooth the noursished profiles
zSmooth(:,1) = smoothdata(nourishedProfile,'gaussian',10,'includenan'); % 50m smoothing
zSmooth(:,2) =  smoothdata(nourishedProfile,'gaussian',4,'includenan'); %20m smoothing
zSmooth(:,3) = nourishedProfile; 

% chunk out areas for smoothing
% let x=-300 to x = -50 have 50 m smoothing
offshore = find(D0045crs>-450 & D0045crs<-20);
zSmooth(offshore,2:3) = NaN;
smoothNourishedProfile = nanmean(zSmooth,2);

% plot(D0045crs,smoothNourishedProfile,'r','linewidth',2)


figureDir = '~/GoogleDriveUCSD/MOPS_Climatology/';
figureName = [figureDir 'IB_CortezProfile.jpeg'];

tt = t(nourished);

text(-650,4,...
    ['Representative Nourished Profile ' datestr(tt(1)) ' - '  datestr(tt(end))],...
    'fontsize',16,'color','red')
% 
% print('-djpeg',figureName,'-r300')

%% save data for SWASH runs
z = smoothNourishedProfile;
x = D0045crs; 


indStart = find(~isnan(z),1);
indStop = find(~isnan(z),1,'last');

z = z(indStart:indStop);
x = x(indStart:indStop);
x = [xEroded(1) x];
z = [zEroded(1); z];

xx = [min(x):2:-202 -200:1:-20 -19:0.5:0 0.25:0.25:max(x)];
zz = interp1(x,z,xx);

bathy.x = xx;
bathy.z = zz;
bathy.dx = 'variable';
bathy.description  = 'MOP line D0045, nourished profile';

plot(x,z,'r','linewidth',2)

% save('~/Documents/Cortez/data/processed/D0045profile_nourished.mat','bathy')
%% Make modern profile using IB lidar scan information
load('~/GoogleDriveUCSD/Cortez/data/lidarscan01_processed.mat','Processed');
plot(Processed.range,Processed.foreshore,'g','linewidth',2)
%% 

h = plot(D0045crs,D0045z(:,85),'b','linewidth',2);
title(datestr(t(85)))

%%
lidarprofile = interp1(Processed.range,Processed.foreshore,D0045crs);
lidarprofile(D0045crs<17) = NaN;

lidarprofile(D0045crs<-50) = smoothErodedProfile(D0045crs<-50);
lidarprofile(D0045crs==0) = D0045z(D0045crs==0,85);
lidarprofile = inpaint_nans(lidarprofile);


zSmoothlidar(:,1) = smoothdata(lidarprofile,'gaussian',4,'includenan'); % 50m smoothing

% plot(D0045crs,lidarprofile,'m','linewidth',2)
plot(D0045crs,zSmoothlidar,'m','linewidth',2)

%% save data for SWASH runs
z = zSmoothlidar;
x = D0045crs; 


indStart = find(x == xEroded(1));
indStop = find(z>10,1);
% 
z = z(indStart:indStop);
x = x(indStart:indStop);
% x = [xEroded(1) x];
% z = [zEroded(1); z];

xx = [min(x):2:-202 -200:1:-20 -19:0.5:0 0.25:0.25:24.75 25:0.1:max(x)];
zz = interp1(x,z,xx);

bathy.x = xx;
bathy.z = zz;
bathy.dx = 'variable';
bathy.description  = 'MOP line D0045, lidar-surveyed eroded profile';

plot(x,z,'m','linewidth',2)

% save('~/Documents/Cortez/data/processed/D0045profile_lidar.mat','bathy')
% save('~/Documents/NeoStockdon/data/D0045profile_lidar.mat','bathy')
%%
% tN = t(indNourishment);
% D0045zN = D0045z(:,indNourishment);
% 
% monthnumN = tN.Month;
% yearnumN = tN.Year;
% 
% mshiftN = nan(1,length(tN));
% % months departure from August, this is beyond klugey but whatever
% for i=1:length(tN)
%     if monthnumN(i)>2
%         mshiftN(i)  = abs(monthnumN(i)-8);
%     elseif monthnumN(i)<=2
%         mshiftN(i) = monthnumN(i)+4;
%     end    
% end


% yearind = [3:9];
% ny = numel(yearind);
% ny = numel(unique(mshift));

% cmap = cmocean('thermal',ny);
% cmap = cmocean('-balance');
% 
% n = 1;
% for i= [7 6 5]
% %     monthind = find(yearnum == years(yearind(n)));
%     monthind = find(mshiftN == i-1);
% 
%     plot(D0045crs,D0045zN(:,monthind),'color',[cmap(round(i/7*256),:) 0.2])
%     hold on
%     n = n+1;
% end
% 
% 
% zNourished = nanmedian(D0045zN(:,ismember(monthnumN,[ 11 12 1 2 3])),2);
% zNourishedSmooth = smoothdata(zNourished,'gaussian',smoothwidth);
% 
% 
% % winter profile
% plot(D0045crs,zNourishedSmooth,'color','k','linewidth',3)
% % summer profile
% % plot(D0045crs,nanmean(D0045zz(:,intersect(monthnum,[6 7 8 9 10])),2),'color',cmap2(2,:),'linewidth',3)
% 





% ax = gca;
% xlims = ax.XLim;
% plot(xlims,MSLnavd88*ones(1,2),':k')
% 
% xlabel('X (m)')
% ylabel('Z (navd88,m)')
% 
% text(xlims(1)+50,MSLnavd88,'MSL','background','w')
% % text(xlims(1)+50,4,['D0045: ' datestr(tt(1)) ' -- '  datestr(tt(end))],'fontsize',16)
% 
% text(xlims(1)+50,4,{'D0045: Median Winter (Nov-Mar) Profiles','with 50m gaussian smoothing'},'color','k','fontsize',16)
% 
% text(xlims(1)+50,3,['NOURISHMENT INFLUENCED: ' datestr(tN(1)) ' - '  datestr(tN(end))],'color','k')
% text(xlims(1)+50,2.5,['NORMALish STATE ' datestr(tt(1)) ' - '  datestr(tt(end))],'color',cmap2(end-1,:))
% 
% 
text(-650,5,...
    ['Representative Lidar+Eroded Profile ' datestr(tt(1)) ' - '  datestr(tt(end))],...
    'fontsize',16,'color','m')

ylim([-11 6])
figureDir = '~/GoogleDriveUCSD/MOPS_Climatology/';
% figureName = [figureDir 'IBprofile_wLidar_zoom.jpeg'];
% % 
% print('-djpeg',figureName,'-r300')

%% alter eroded profile to match lidar profile with 15m depth.
clear all
% Mean Seal Level in NAVD88
NAVD88mllw = 0.058;
MSLmllw = 0.832;
MSLnavd88 = MSLmllw-NAVD88mllw;
MHWmllw = 1.402;
MHWnavd88 = MHWmllw-NAVD88mllw;
eroded = load('~/Documents/NeoStockdon/data/D0045profile_eroded_15.mat');
plot(eroded.bathy.x,eroded.bathy.z)
hold on


%
load('~/GoogleDriveUCSD/Cortez/data/lidarscan01_processed.mat','Processed');
plot(Processed.range,Processed.foreshore,'g','linewidth',2)
% 
% % make x^2 function to raise slope
% T = ones(length(bathy.x),1);
% inflectionpt = find(bathy.x>35,1);
% 
% y = bathy.z;
% y(bathy.x>0) = 1e-3.*(bathy.x(bathy.x>0)).^2.2;
% 
% y(bathy.x>-5 & bathy.x<=0) = 1e-4.*(bathy.x(bathy.x>-5 & bathy.x<=0)).^2;
% 
% xx = bathy.x;
% zz = bathy.z;
% 
% plot(xx,y)
lidar = load('~/Documents/NeoStockdon/data/D0045profile_lidar.mat');
plot(lidar.bathy.x,lidar.bathy.z)
x = eroded.bathy.x;
z = interp1(lidar.bathy.x,lidar.bathy.z,x);

z(x<=-650) = eroded.bathy.z(x<=-650);

% extend to z=20
ind = find(~isnan(z),1,'last');
mdl = fitlm(x(ind-1:ind),z(ind-1:ind));
mb = mdl.Coefficients.Estimate;

newslope = mb(1)+x*mb(2);
xind = knnsearch(newslope(:),20+MSLnavd88);

z(ind:xind) = newslope(ind:xind);

stopInd = find(isnan(z));
z(stopInd) = [];
x(stopInd) = [];

plot(x,z)

bathy.x = x;
bathy.z = z;
bathy.description = 'MOP line D0045, lidar-surveyed eroded profile, 15m depth';
bathy.dx = lidar.bathy.dx;

save('~/Documents/NeoStockdon/data/D0045profile_lidar_15.mat','bathy')



