% clear all

filedir = '/Volumes/group/LiDAR/Waves/20190118_MOP45_CortezRunup/IB_2019_01_17_and_18_Rnup/2019_01_18/Line_scanning/UTM_GMT/';
filename = 'LineScan1_utm11n_navd88_intensity_gmt.mat';

load([filedir filename])
%%
parosX = 487545.49;
parosY = 3603947.41;
%% take off base x and base y for easier units
% for ifile = 3:5

x = xyzit(:,1);
y = xyzit(:,2);
z = xyzit(:,3);
amp = xyzit(:,4);
t = xyzit(:,5);
%%
x0 = 4.8757355481e5;
y0 = 3.60394954046e6 ;

x = x-x0;
y = y-y0;
r=sqrt(x.^2+y.^2); %range

parosx = (parosX-x0);
parosy = (parosY-y0);

% put paros on lidar line assuming lidar direction is shorenormal and beach
% is alongshore uniform.

ploc = [parosx parosy];
lidarloc = [x y];
lidarloc = lidarloc(1:10000:end,:);

for i=1:length(lidarloc)
X = [lidarloc(i,:); ploc];
parosmin(i) = pdist(X);
end

indmin = find(parosmin==min(parosmin));

parosProjection = lidarloc(indmin,:);

parosR = norm(parosProjection);
parosZ = 1.56;
%%
clf
plot(r,'.')
%%
% for ifile = 3:5
% xgrid = [20:0.1:120];

% for ifile = 9
% xgrid = [0:0.1:150];
% 
% ygrid = [780:0.1:820];

%% loop on t
scanN = find(abs(diff(t)>0.8e-6));
%%
% scanN = find(abs(diff(scangle))>200);

for i=1:length(scanN)
tt(i) = t(scanN(i));
end
%%
indstart = nan(length(scanN)-1,1);
%find when to start the scan
for i=1:length(scanN)-1
    sctemp = y(scanN(i):scanN(i+1));
%     ltemp = min(sctemp);
    minloc = find(sctemp == max(sctemp),1);
    indstart(i) = scanN(i)+minloc-1;
end

% add first scan
sctemp = y(1:scanN(1));
minloc = find(sctemp == max(sctemp),1);
indstart = [scanN(1)+minloc-1; indstart];

%%
clf
plot(y,'.')
hold on
plot(scanN,y(scanN),'or')
plot(indstart,y(indstart),'og')
% amp = get_intensity(c);
%%

tot_scans = numel(scanN);
scanDIFF =scanN-indstart;
ind_scanstart = indstart;
numPTSscan = scanDIFF;
I = tot_scans -2;
J = max(scanDIFF);

Tmat = NaN(I,J);
Rmat = NaN(I,J);
Zmat = NaN(I,J);
Amat = NaN(I,J);

for i=1:I
    Tmat(i,1:scanDIFF(i)) = t(ind_scanstart(i):(ind_scanstart(i)+numPTSscan(i)-1));
    Rmat(i,1:scanDIFF(i)) = r(ind_scanstart(i):(ind_scanstart(i)+numPTSscan(i)-1));
    Zmat(i,1:scanDIFF(i)) = z(ind_scanstart(i):(ind_scanstart(i)+numPTSscan(i)-1));
    Amat(i,1:scanDIFF(i)) = amp(ind_scanstart(i):(ind_scanstart(i)+numPTSscan(i)-1));

    
    ind_good = find(~isnan(Tmat(i,:)));
    Tmat(i,1:length(ind_good)) = fliplr(Tmat(i,ind_good));
    Tmat(i,1+length(ind_good):J) = NaN;
    Amat(i,1:length(ind_good)) = fliplr(Amat(i,ind_good));
    Amat(i,1+length(ind_good):J) = NaN;
    Rmat(i,1:length(ind_good)) = fliplr(Rmat(i,ind_good));
    Rmat(i,1+length(ind_good):J) = NaN;
    Zmat(i,1:length(ind_good)) = fliplr(Zmat(i,ind_good));
    Zmat(i,1+length(ind_good):J) = NaN;
end

%%
[I,J] = size(Rmat);
xi_interp = [0:0.1:50]; % HARD CODE for region we care about

Zinterp = nan(I,length(xi_interp));
for m=3:I
    X = reshape(Rmat(m,:),1,J);
    Z = reshape(Zmat(m,:),1,J);
    [C,ia,ic] = unique(X,'stable');
    X = X(ia);
    Z = Z(ia);
    if sum(~isnan(X))>2
    Zinterp(m,:) = interp1(X(~isnan(X)),Z(~isnan(X)),xi_interp);
    end
end

%%
[M,~]=size(Rmat);
xi = xi_interp;
dxi=mean(diff(xi))/2;
Zinterp2=nan(M,length(xi));
Ainterp=nan(M,length(xi));
%%
for a=1:M
    r=Rmat(a,:);
    z=Zmat(a,:);
    amp=Amat(a,:);
    for i=1:length(xi)
        ind=r>=xi(i)-dxi & r<=xi(i)+dxi;
        if ~isempty(z(ind))
            Zinterp2(a,i)=median(z(ind));
            Ainterp(a,i)=median(amp(ind));
        end
    end
end

xi_interp = -xi_interp;

%denoise
%% save data
Processed.range = xi_interp;
Processed.Zinterp = Zinterp;
Processed.Zinterp2 = Zinterp2;
Processed.Ainterp = Ainterp;
Processed.t = Tmat(:,1)';


% save('../mat/ProcessedLidar_5','filename','Processed','xyzti');
% save('../mat/ProcessedLidar_5','filename','xyzti');

%% 
clf
pcolor(Processed.range,Processed.t,Processed.Zinterp2);
shading flat
xlabel('m from lidar')
datetick('y','MM:SS')
ylabel('Time (MM:SS)')
title(filename, 'Interpreter', 'none');
hc = colorbar;
caxis([-0.5 5])
hc.Label.String = 'Z (m, NAVD88)';
%%
hFig = gcf;
hFig.PaperUnits = 'inches';
hFig.PaperSize = [5 7];
hFig.PaperPosition = [0 0 5 7];
% print(hFig, '-djpeg', [savedir filename '.jpg'],'-r300');

% caxis([0.5 3])
% %% 
% hold on
% load sensorlocs.mat
% 
% plot([Paros(1).x Paros(1).x],[Tmat(2,1) Tmat(end,1)],'r')
% plot([Paros(2).x Paros(2).x],[Tmat(2,1) Tmat(end,1)],'r')
% plot([Paros(3).x Paros(3).x],[Tmat(2,1) Tmat(end,1)],'r')
% plot([Paros(4).x Paros(4).x],[Tmat(2,1) Tmat(end,1)],'r')
% plot([Paros(5).x Paros(5).x],[Tmat(2,1) Tmat(end,1)],'r')
% 
% 
%
% %%
% hFig = figure
% hFig.PaperSize = [7 3];
% hFig.PaperPosition = [0 0 7 3];
%%
indvect = find(vect.time(1,:) >Tmat(2,1),1);
clf
ax(1) = subplot(3,1,1:2);
plot(xi_interp(1:351),nanmin(Zinterp2(:,1:351)),'r')
hold on
scatter(-parosR,parosZ,'+m','linewidth',2)
plot([-50 0],[1 1]*mean(vect.press(:,indvect))-6.1,'b')
n = 1;
start = 3000;
stop = 3600;
xlabel('m from lidar')
ylabel('Z (m, NAVD88)')
    
    
ax(2) = subplot(3,1,3);
indx = find(xi_interp < parosx,1);
indt = find(paros.t(1,:) >Tmat(2,1),1);

tdate = datetime(Tmat(:,1),'ConvertFrom','datenum');

plot(tdate,Zinterp(:,indx),'color',[0.3 0.3 0.3])
hold on
plot(paros.date(:,indt),paros.P(:,indt)+1.599,'m')
plot(tdate,Zinterp2(:,indx),'.k')
xlim([tdate(start) tdate(stop)])
ylabel('Z (m,navd88)')
ylim([0.75 3])

for i = start:stop
    axes(ax(1))
    h(4:6) = plot(xi_interp,Zinterp2(i-3:i-1,:),'.','color',[0.7 0.7 0.7]);
    h(2) = plot(xi_interp,Zinterp2(i,:),'.k');
    h(3) = text(-45,5.5,[datestr(Processed.t(i)) ' UTC']);
    xlim([-50 0])
    ylim([-.5 6])
    pause(0.1)
    
    
    axes(ax(2))
    th = plot([tdate(i) tdate(i)],[0 3],'g');
    
    drawnow
    F = getframe(gcf);
    [im{n},map] = frame2im(F);
    if i<stop
        axes(ax(1))
        delete(h)
        axes(ax(2))
        delete(th)
    end
    n= n+1;
end



% 
% 
% %%
% currdir = pwd;
cd /Users/juliafiedler/Desktop
gifname = [filename(1:end-4) '.gif'];
for idx = 1:length(im)
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,gifname,'gif','LoopCount',Inf,'DelayTime',0);
    else
        imwrite(A,map,gifname,'gif','WriteMode','append','DelayTime',0);
    end
end
% cd(currdir)


% for i=1:length(tt)-1
% xtemp = x(scanN(i)+1:scanN(i+1));
% xu = find(unique(xtemp));
% ztemp = z(scanN(i)+1:scanN(i+1));
% xtemp = xtemp(xu);
% ztemp = ztemp(xu);
% 
% zinterp(i,:) = interp1(xtemp,ztemp,xgrid,'nearest');
% end