% this really needs to be spruced up.
% made into a function or something
% process_runupStatsMOPS(MOPname,yind)
clear
MOPname = 'D0045';
MOPdata = load('~/Documents/NeoStockdon/data/binnedMOPS.mat','ZBIN','MOP','Ho','LoC');


for yind=47
    inputfiledir = [ '/Volumes/jfiedler/NeoStockdon/D0045_Eig015/' num2str(yind,'%02.0f')  ];
    outputfiledir = [ '~/Documents/NeoStockdon/processed/D0045_Eig015/' num2str(yind,'%02.0f')  ];
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
    
    % run the loop backwards to preallocate by default
    for n=nfol:-1:1
        n
        modDir{n} = [filelist(n).folder '/' filelist(n).name '/'];
%         if ismember(yind,33:35)
%             modDir{n} = [filelist(n).folder '/'];
%         end
        specname = dir([modDir{n} '2*.spec']); %<-- this will not work with 2D!! TODO: FIX THIS
        
        date1 = datetime([specname.name(1:11) '0000']);
        tind(n) = find(MOPdata.MOP.time == date1);
        options(n).Tm = 1/MOPdata.MOP.fp(tind(n));
        
        %% get water level (tide) from sws input file
        
        fname = dir([modDir{n} '*.sws']);
        c = fileread([modDir{n} fname(end).name]);
        t = regexp(c,'wlev = \w*[.]\d{3}','start');
        options(n).wlev = str2double(c(t+7:t+11));
        options(n).windowlength = 7;
        
        try
            [Spec(n),Info(n),Bulk(n),Tseries(n),Sur(n)] = get_runupStats(modDir{n},date1,options(n));
            Sig(n) = Bulk(n).swashparams(1);
            Sinc(n) = Bulk(n).swashparams(2);
            etaRunup(n) = Bulk(n).swashparams(3);
            beta(n) = Bulk(n).beta;
            Eig(n) = Bulk(n).Eig;
        catch ME % TODO: FIx this! THis is the worst.
            warning('No runup stats, assigning nans to Bulk stats')
            Sig(n) = nan;
            Sinc(n) =nan;
            etaRunup(n) = nan;
            beta(n) = nan;
            Eig(n) = nan;
        end
    end
    %%
    clear c t fname date1 n
    filestem = '/Users/juliafiedler/Documents/NeoStockdon/processed/D0045_Eig015/';
    %         fname = [filestem num2str(yind,'%02.0f') '/bulkrunup_fspread_' num2str(xind,'%02.0f') '.mat'];
    fname = [filestem num2str(yind,'%02.0f') '/bulkrunup.mat'];
    
    if exist('Bulk','var')
        save(fname, 'Bulk','Info','MOPname','Spec','Sur','Tseries','filelist','modDir','options','specname','yind','tind')
    else
        save(fname,'ME','filelist','modDir','options','specname','yind','tind')
    end
    
    
    Ho = MOPdata.Ho(tind)';
    Lo = MOPdata.LoC(tind)'; % <---***ATTENTION THIS IS THE CENTROID WAVELENGTH****
    
    S = sqrt(Sig.^2+Sinc.^2);
    R2_stat = 1.1*(etaRunup + S/2);
    
    
    [R2_S06,Sig_S06,Sinc_S06,etaRunup_S06,sqrtHoLo] = get_Stockdon2006params(Ho,Lo,beta);
    
    
    fname = [filestem num2str(yind,'%02.0f') '/runupstats.mat'];
    save(fname,'R2_stat','R2_S06','S','Sig','Sig_S06','Sinc','Sinc_S06','etaRunup_S06','etaRunup','beta','Ho','Lo','sqrtHoLo','modDir','tind','yind','Eig')
    
    clearvars -except MOPdata yind MOPname
    
    
end



