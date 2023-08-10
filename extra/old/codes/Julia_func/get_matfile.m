function [res,datahour,botlev,wlev,Wlev,X,M,N, bulkinfo] = get_matfile(folderdir,filedir,casenum)



filename = 'agate.mat';


[res,~] = swash_loadMatFile([folderdir filedir filename]);
filename = sprintf('agate%02.0f.sws',casenum);
fid = fopen([folderdir filedir filename]);
for i=1:5
    fgetl(fid);
end
datahour = datenum(fgetl(fid),'$ DD-mmm-yyyy HH:MM');
bulkinfo = fgetl(fid); bulkinfo = bulkinfo(3:end);
wlev = fgetl(fid); wlev = wlev(10:end-1);
wlev = str2double(wlev);

% Wlev = squeeze(res.Watlev(1,:,end-7199:end));
Wlev = squeeze(res.Watlev(1,:,:));

[M,N] = size(Wlev);
X = res.Xp(1,:); 

botlev = -(res.Botlev-wlev);

clear i fid ans