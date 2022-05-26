function get_model_sectormoments(folder,savename,casenumber,directionality,inorm, flim)

% THIS USES NORM TO TOTAL VARIANCE FOR PLOTTING PURPOSES
if ~exist('inorm','var')
  inorm = 0;
end

load ~/SWASH/data/PUV_for_Physics Agate
for i=1:19
    t(i) = Agate(i).t;
end

clearvars -except t savename folder casenumber directionality inorm flim

t(20) = datenum(2013,10,1,19,0,0);
t(21) = datenum(2013,10,1,22,0,0);
t(22) = datenum(2013,10,21,13,0,0);

%%

if nargin<3
    casen = 1;
    while casen<22
        if casen==4
            casen = 6;
        end
        datahour = t(casen);
        filedir = sprintf('agate_%02.0f/',casen);
        folderdir = [folder filedir];
        %     folderdir = ['~/Desktop/agate_09_weak/']
        filename= 'agate.mat';
        
        [res,vars] = swash_loadMatFile([folderdir filename]);
        
        xarray = res.Xp(1,:);
        bot = -res.Botlev(1,:);
        wlev = squeeze(res.Watlev(1,:,:));
        wlev = wlev-bot'*ones(1,length(wlev));
        wlev = wlev';
        U = squeeze(res.vel_x(1,:,:));
        %%
        T = 7100/2;
        sig = wlev(1:7100,:);
        
        sigIN = get_directionality(sig,U,directionality);
        
        %%
%             flim = [0.004 0.025 0.04 0.25];
%         flim = [0.004 0.04 0.25];
        [m2,m3] = sectormoments( T, sigIN , flim,inorm);
        %%
        biphase.self = atan2(imag(m3.self),real(m3.self));
        biphase.sum = atan2(imag(m3.sum),real(m3.sum));
        biphase.dif = atan2(imag(m3.dif),real(m3.dif));
        %%
        
        biphase_self(casen,:,:)  = biphase.self;
        biphase_sum(casen,:,:)  = biphase.sum;
        biphase_dif(casen,:,:)  = biphase.dif;
        %
        skew_self(casen,:,:) = real(m3.self);
        skew_sum(casen,:,:) = real(m3.sum);
        skew_dif(casen,:,:) = real(m3.dif);
        skew_tot(casen,:,:) = real(m3.tot);
        %
        asym_self(casen,:,:) = imag(m3.self);
        asym_sum(casen,:,:) = imag(m3.sum);
        asym_dif(casen,:,:) = imag(m3.dif);
        
        bicoh_self(casen,:,:) = sqrt(dot(m3.self,m3.self,3));
        bicoh_sum(casen,:,:) = sqrt(dot(m3.sum,m3.sum,3));
        bicoh_dif(casen,:,:) = sqrt(dot(m3.dif,m3.dif,3));
        
        casen = casen+1;
    end
end

if nargin>2
    casen = casenumber;
    datahour = t(casen);
    %     filedir = sprintf('agate_%02.0f/',casen);
    %     folderdir = [folder filedir];
    %     folderdir = ['~/Desktop/agate_09_weak/']
    filename= 'agate.mat';
    
    [res,vars] = swash_loadMatFile([folder filename]);
    
    xarray = res.Xp(1,:);
    bot = -res.Botlev(1,:);
    wlev = squeeze(res.Watlev(1,:,:));
    wlev = wlev-bot'*ones(1,length(wlev));
    wlev = wlev';
    
    %%
    [M,N] = size(wlev);
    u_vel = nan(10,N,M);
    for i=1:10
        varname = sprintf('res.vel_k%02.0f_x',i);
        eval(['u_vel(i,:,:)  =' varname ';']);
    end
    %%
    U = mean(u_vel);
    U = squeeze(U(1,:,:));
    %%
    T = 7100/2;
    sig = wlev(1:7100,:);
    
    sigIN = get_directionality(sig,U,directionality);
    %%
%         flim = [0.004 0.025 0.04 0.25];
%     flim = [0.004 0.04 0.25];
    [m2,m3] = sectormoments( T, sigIN , flim,inorm);
    %%
    biphase.self = atan2(imag(m3.self),real(m3.self));
    biphase.sum = atan2(imag(m3.sum),real(m3.sum));
    biphase.dif = atan2(imag(m3.dif),real(m3.dif));
    %%
    
    biphase_self  = biphase.self;
    biphase_sum  = biphase.sum;
    biphase_dif  = biphase.dif;
    %
    skew_self = real(m3.self);
    skew_sum = real(m3.sum);
    skew_dif = real(m3.dif);
    skew_tot = real(m3.tot);
    %
    asym_self = imag(m3.self);
    asym_sum = imag(m3.sum);
    asym_dif = imag(m3.dif);
    
    
    bicoh_self = sqrt(dot(m3.self,m3.self,3));
    bicoh_sum = sqrt(dot(m3.sum,m3.sum,3));
    bicoh_dif = sqrt(dot(m3.dif,m3.dif,3));
    
    
end

    function sigIN = get_directionality(sig,U,directionality)
        switch directionality
            case 'shoreward'
                disp('Removing reflection')
                %% Remove reflection
                H = mean(sig)';
                u = U(:,1:7100)';
                fmin = 0.001; fmax=0.04;
                for i=1:length(bot)
                    [sigIN(:,i),sigOUT(:,i)] = RemoveReflection( sig(:,i), u(:,i), H(i),  0.5,'method','sher','flim', [fmin,fmax],'g',9.81,'trend',true,'depav',true);
                end
            case 'all'
                sigIN = sig;
                disp('Not removing reflection')
            otherwise
                warning('Directionality not correctly specified.')
        end
    end


clearvars -except t skew* biphase* bicoh* asym* savename xarray flim
clear biphase
save(savename)
disp(['Saved to ' savename ' : ' datestr(now)])
end