% this really needs to be spruced up.
% made into a function or something
% process_runupStatsMOPS(MOPname,yind)
clear
addpath  /Users/juliafiedler/Documents/SWASH_runup/mfiles/functions
addpath  /Users/juliafiedler/Documents/SWASH_runup/mfiles/processdata/


MOPname = 'D0045';
MOPdata = load('~/Documents/NeoStockdon/data/binnedMOPS.mat','ZBIN','MOP','Ho','LoC');



% inputfiledir = '~/Documents/NeoStockdon/models/D0045_Eig015/';
inputfiledir = '~/Documents/NeoStockdon/models/Eig015_lidar/';

outputfiledir = '~/Documents/NeoStockdon/processed/Eig015_lidar/';
thresholdheight = 0.03;
bulkrunupname = 'bulkrunup_3cm.mat';
runupstatsname = 'runupstats_3cm.mat';


if ~exist(outputfiledir,'dir')
    mkdir(outputfiledir)
end
%TODO: loop through folders in models to get different input dates and make
%new structure with all runupStats
%         filelist = dir([inputfiledir '/' num2str(xind,'%02.0f') '/2*']);
%     if yind>=32
%         filelist = dir([inputfiledir '/*/2*']);
%     elseif yind<32
%         filelist = dir([inputfiledir '/2*']);
%     end
filelist = dir([inputfiledir '/2*']);

nfol = length(filelist);
%%
% run the loop backwards to preallocate by default
for n=nfol
    n
    modDir{n} = [filelist(n).folder '/' filelist(n).name '/'];
    %         if ismember(yind,33:35)
    %             modDir{n} = [filelist(n).folder '/'];
    %         end
    specname = dir([modDir{n} '2*.spec']); %<-- this will not work with 2D!! TODO: FIX THIS
    
    date1 = datetime([specname.name(1:11) '0000']);
    tind(n) = find(MOPdata.MOP.time == date1);
    options(n).Tm = 1/MOPdata.MOP.fp(tind(n));
    options(n).threshold = thresholdheight;
    options(n).totaltime = 60;

    
    %% get water level (tide) from sws input file
    
    fname = dir([modDir{n} '*.sws']);
    c = fileread([modDir{n} fname(end).name]);
    t = regexp(c,'wlev = \w*[.]\d{3}','start');
    options(n).wlev = str2double(c(t+7:t+11));
    options(n).windowlength = 10;
    
%     try
        [Spec(n),Info(n),Bulk(n),Tseries(n),Sur(n)] = get_runupStats(modDir{n},date1,options(n));
        Sig(n) = Bulk(n).swashparams(1);
        Sinc(n) = Bulk(n).swashparams(2);
        etaRunup(n) = Bulk(n).swashparams(3);
        beta(n) = Bulk(n).beta;
        Eig(n) = Bulk(n).Eig;
        Ess(n) = Bulk(n).Ess;
%     catch ME % TODO: FIx this! THis is the worst.
%         warning('No runup stats, assigning nans to Bulk stats')
%         disp(modDir{n})
%         Sig(n) = nan;
%         Sinc(n) =nan;
%         etaRunup(n) = nan;
%         beta(n) = nan;
%         Eig(n) = nan;
%         Ess(n) = nan;
%     end
end

% run the loop backwards to preallocate by default
parfor n=1:nfol-1
    n
    modDir{n} = [filelist(n).folder '/' filelist(n).name '/'];
    %         if ismember(yind,33:35)
    %             modDir{n} = [filelist(n).folder '/'];
    %         end
    specname = dir([modDir{n} '2*.spec']); %<-- this will not work with 2D!! TODO: FIX THIS
    
    date1 = datetime([specname.name(1:11) '0000']);
    tind(n) = find(MOPdata.MOP.time == date1);
    options(n).Tm = 1/MOPdata.MOP.fp(tind(n));
    options(n).threshold = thresholdheight;
    options(n).totaltime = 60;

    
    %% get water level (tide) from sws input file
    
    fname = dir([modDir{n} '*.sws']);
    c = fileread([modDir{n} fname(end).name]);
    t = regexp(c,'wlev = \w*[.]\d{3}','start');
    options(n).wlev = str2double(c(t+7:t+11));
    options(n).windowlength = 10;
    
    dateN(n) = date1;
end
 %%
 disp('moving on')

parfor n=1:nfol-1
    n
    try
        [Spec(n),Info(n),Bulk(n),Tseries(n),Sur(n)] = get_runupStats(modDir{n},dateN(n),options(n));
        Sig(n) = Bulk(n).swashparams(1);
        Sinc(n) = Bulk(n).swashparams(2);
        etaRunup(n) = Bulk(n).swashparams(3);
        beta(n) = Bulk(n).beta;
        Eig(n) = Bulk(n).Eig;
        Ess(n) = Bulk(n).Ess;
    catch ME % TODO: FIx this! THis is the worst.
        warning('No runup stats, assigning nans to Bulk stats')
        disp(modDir{n})
        Sig(n) = nan;
        Sinc(n) =nan;
        etaRunup(n) = nan;
        beta(n) = nan;
        Eig(n) = nan;
        Ess(n) = nan;
    end
end
%%
clear c t fname date1 n
% filestem = '/Users/juliafiedler/Documents/NeoStockdon/processed/D0045_Eig015/';
%         fname = [filestem num2str(yind,'%02.0f') '/bulkrunup_fspread_' num2str(xind,'%02.0f') '.mat'];
fname = [outputfiledir bulkrunupname];

if exist('Bulk','var')
    save(fname, 'Bulk','Info','MOPname','Spec','Sur','Tseries','filelist','modDir','options','specname','tind')
else
    save(fname,'ME','filelist','modDir','options','specname','tind')
end

%% FIX THIS SECTION!! TODO
Ho = MOPdata.Ho(tind)';
Lo = MOPdata.LoC(tind)'; % <---***ATTENTION THIS IS THE CENTROID WAVELENGTH****

S = sqrt(Sig.^2+Sinc.^2);
R2_stat = 1.1*(etaRunup + S/2);


[R2_S06,Sig_S06,Sinc_S06,etaRunup_S06,sqrtHoLo] = get_Stockdon2006params(Ho,Lo,beta);


fname = [outputfiledir runupstatsname];
save(fname,'R2_stat','R2_S06','S','Sig','Sig_S06','Sinc','Sinc_S06','etaRunup_S06','etaRunup','beta','Ho','Lo','sqrtHoLo','modDir','tind','Eig','Ess')

clearvars -except MOPdata yind MOPname





