function get_H_bathytests(inFileDir,outFileDir,datahour)
%% get all swash wave heights across surfzone

% volumedir = '/Volumes/jfiedler/SWASH_runup';

% inFileDir = [volumedir '/mat/bathytest0_005f/'];

% outFileDir = [volumedir 'mat/processed/'];
outFileName = 'H_bathytests.mat';

% datahour = datenum(2013,10,21,13,0,0);
windowlength = 5; % 5 minutes window length

% TODO: possibly automate this if we're ever going to change the Hz
% collection.
nfft = windowlength*60*5; % windowlength(minutes) * 60 seconds * 5 Hz
filelist = dir([inFileDir '*mat']);

for i=1:length(filelist)
    disp(['Running ' filelist(i).name])
    try
    S(i) = getHcase([inFileDir num2str(i,'%02.0f') '.mat'],nfft,'dfdx');
    S(i).Hinc_S(5)
    catch
    end
end
    

function casename = getHcase(folderdir,nfft,test)
[casename.Hinc_S,casename.Hig_S,~,~,~,casename.X_S, casename.datahour, ~,~,casename.CI_Hinc_S,casename.CI_Hig_S]...
    = get_H_spatial(folderdir,22,'total','ebalance',nfft,test);
%     get_H_spatial(folderdir,casenum,method,type,nfft,test)
casename.casename = folderdir;
end

%% get all paros/PUV wave heights across surfzone
[P.X, P.Hinc, P.Hig, P.nameSensors] = get_H_paros(datahour);
P.datahour = datetime(datahour,'ConvertFrom','datenum');

save([outFileDir outFileName],'S','P')

end