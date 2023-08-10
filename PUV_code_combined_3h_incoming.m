%% PUV processing - from raw data files to all stats, boundwave and bispectra computed
clear all
tic
% input directory where data is - one deloyment per folder
local_dir = uigetdir(pwd,'Local Directory'); % Volumes/LANGE_Passport/PUV
cd(local_dir)

code_dir = fullfile(local_dir, 'CODES');
data_dir = fullfile(local_dir, 'DATA', 'PUV');
addpath(genpath(code_dir))
addpath(genpath(data_dir))
mkdir([data_dir '/Level1_QC/'])
mkdir([data_dir '/Level2_QC/'])
mkdir([data_dir '/Level3_QC/'])
 
%% Pull deploy_notes table from Excel and sort by MOP number
% read in table of filename that contains 'Notes'
files = dir(data_dir);files([files.isdir]==1) = [];
for ii = 1:length(files);names{ii}=files(ii).name;end
deploy_notes = readtable(string(names(contains(names, 'Notes'))));

[~,ii]=sort(deploy_notes.MOP);
deploy_notes=deploy_notes(ii,:);
aa=tabulate(deploy_notes.MOP); 

% sort by MOP number
if length(find(aa(:,2)>1))>0
    mm = find(aa(:,2)>1); % if multiple sensors on single mop line
    for jj = 1:length(mm)
        ll = mm(jj);
        nn = find(aa(ll,1)==deploy_notes.MOP);
        [~,ii]=sort(deploy_notes.Depth(nn));
        deploy_notes(nn,:)=deploy_notes(nn(ii),:);
    end
end
% 
% % remove any . files or spare things
% for nn = length(files):-1:1 
%     if files(nn).name(1) == '.' | files(nn).isdir == 0 | contains(files(nn).name, 'Level')
%         files(nn)=[];
%     end
% end
% %% remove LosPen stuff for now
% for nn = size(deploy_notes,1):-1:1
%     if contains(char(string(deploy_notes{nn,3})), 'LosPen') | contains(char(string(deploy_notes{nn,3})), 'Del Mar') | contains(char(string(deploy_notes{nn,3})), 'IB') | deploy_notes{nn,4}==15
%         deploy_notes(nn,:)=[];
%     end
% end

%% Raw Process Data - Level 1
for ii = 1:size(deploy_notes,1)
    directory = [data_dir '/Raw_DATA']
    
    filename = [char(string(deploy_notes.Location(ii))) '_' char(string(deploy_notes.Date(ii))) '_MOP' num2str(deploy_notes.MOP(ii)) '_' num2str(deploy_notes.Depth(ii)) 'm']
    % input lat/lon of sensor
    LATLON = [deploy_notes.Latitude(ii) deploy_notes.Longitude(ii)];
    % input compass heading of beam 1 of ADV
    rot_angle = [deploy_notes.Heading(ii)];
    % input clock drift between beginning and end of deployment
    clockdrift = [deploy_notes.ClockDrift(ii)];
    [PUV] = PUV_raw_process(directory, filename, LATLON, rot_angle, clockdrift)
    save([data_dir '/Level1_QC/' filename '_processed.mat'], 'PUV', 'filename')
    
end

%% Compute PUV_process - Level 2
for ll = 1:size(deploy_notes,1)
    tic
    clearvars -except files deploy_notes names ll data_dir code_dir
    disp('New run')
    %directory = ['MOP' num2str(deploy_notes.MOP(ll)) '-' num2str(deploy_notes.Depth(ll)) 'm'];
    filename = [char(string(deploy_notes.Location(ll))) '_' char(string(deploy_notes.Date(ll))) '_MOP' num2str(deploy_notes.MOP(ll)) '_' num2str(deploy_notes.Depth(ll)) 'm']
    directory = [data_dir '/Level1_QC/']
    load([directory '/' filename '_processed.mat'], 'PUV')
    aa = char(string(split(filename, '_')));
    PUV.mop = aa(3,4:6);
    doffp = deploy_notes.Doffp(ll)/100;

    %% Rotate to shorenormal coordinates (MOP)
    [PUV.ShorenormalCoord.U, PUV.ShorenormalCoord.V, shorenormal] = rotate_shorenormal(PUV.BuoyCoord.U, PUV.BuoyCoord.V, PUV.mop);
    
    % 66cm sand to top of pressure port on deployment. 78cm on recovery. 
    % doffp = 0.76; % splitti   ng the difference between the two
    % take doffp from Deployment Notes
    doffu = doffp+0.371; %adding 37.1 cm to u,v sample point
    PUV_work = PUV;
    %% Interp NaNs if less than 2sec gap
    disp('interp 2sec')
    PUV_work.P = fillmissing(PUV.P,'linear','SamplePoints',PUV.time,'MaxGap',seconds(2));
    PUV_work.BuoyCoord.U = fillmissing(PUV.BuoyCoord.U,'linear','SamplePoints',PUV.time,'MaxGap',seconds(2));
    PUV_work.BuoyCoord.V = fillmissing(PUV.BuoyCoord.V,'linear','SamplePoints',PUV.time,'MaxGap',seconds(2));
    PUV_work.BuoyCoord.W = fillmissing(PUV.BuoyCoord.W,'linear','SamplePoints',PUV.time,'MaxGap',seconds(2));
    PUV_work.ShorenormalCoord.U = fillmissing(PUV.ShorenormalCoord.U,'linear','SamplePoints',PUV.time,'MaxGap',seconds(2));
    PUV_work.ShorenormalCoord.V = fillmissing(PUV.ShorenormalCoord.V,'linear','SamplePoints',PUV.time,'MaxGap',seconds(2));
    PUV_work.T = fillmissing(PUV.T,'linear','SamplePoints',PUV.time,'MaxGap',seconds(2));
    %% Check for missing or too many data points
    % check that only 7200 values recorded every hour - add NaN if a value
    % is skipped. Remove value if extra value recorded
    disp('check for missing points')
    it=find(floor(second(PUV_work.time))==0 & minute(PUV_work.time) == 0);
    for ii = 2:length(it)
        if diff(it(ii-1:ii)) == 1
            it(ii)=NaN;
        end
    end
    it(isnan(it))=[];

    for ii = 1:length(it)-1
        if it(ii+1)-it(ii) == 7201 % extra value
            PUV_work.time(it(ii+1)-1)=[];
            PUV_work.P(it(ii+1)-1)=[];
            PUV_work.BuoyCoord.U(it(ii+1)-1)=[];
            PUV_work.BuoyCoord.V(it(ii+1)-1)=[];
            PUV_work.BuoyCoord.W(it(ii+1)-1)=[];
            PUV_work.ShorenormalCoord.U(it(ii+1)-1)=[];
            PUV_work.ShorenormalCoord.V(it(ii+1)-1)=[];
            it(ii+1:end) = it(ii+1:end)-1;
           
        elseif it(ii+1)-it(ii) == 7199 % missing a value
            PUV_work.time = [PUV_work.time(1:it(ii)+7200-1); PUV_work.time(it(ii)+7200-1); PUV_work.time(it(ii)+7200:end)];
            PUV_work.P = [PUV_work.P(1:it(ii)+7200-1); NaN; PUV_work.P(it(ii)+7200:end)];
            PUV_work.BuoyCoord.U = [PUV_work.BuoyCoord.U(1:it(ii)+7200-1); NaN; PUV_work.BuoyCoord.U(it(ii)+7200:end)];
            PUV_work.BuoyCoord.V = [PUV_work.BuoyCoord.V(1:it(ii)+7200-1); NaN; PUV_work.BuoyCoord.V(it(ii)+7200:end)];
            PUV_work.BuoyCoord.W = [PUV_work.BuoyCoord.W(1:it(ii)+7200-1); NaN; PUV_work.BuoyCoord.W(it(ii)+7200:end)];
            PUV_work.ShorenormalCoord.U = [PUV_work.ShorenormalCoord.U(1:it(ii)+7200-1); NaN; PUV_work.ShorenormalCoord.U(it(ii)+7200:end)];
            PUV_work.ShorenormalCoord.V = [PUV_work.ShorenormalCoord.V(1:it(ii)+7200-1); NaN; PUV_work.ShorenormalCoord.V(it(ii)+7200:end)];
            it(ii+1:end) = it(ii+1:end)+1;
           
        end
    end
    %% Split data into 3 hour segments
    disp('split into 3h segments')
    hourlen = (60*60*PUV_work.fs)*3;
    segtotal = floor(length(PUV_work.time)/hourlen);
    PUV_work.time(hourlen*segtotal)

    PUV_work.time(segtotal*hourlen+1:end)=[];
    PUV_work.P(segtotal*hourlen+1:end)=[];
    PUV_work.BuoyCoord.U(segtotal*hourlen+1:end)=[];
    PUV_work.BuoyCoord.V(segtotal*hourlen+1:end)=[];
    PUV_work.BuoyCoord.W(segtotal*hourlen+1:end)=[];
    PUV_work.ShorenormalCoord.U(segtotal*hourlen+1:end)=[];
    PUV_work.ShorenormalCoord.V(segtotal*hourlen+1:end)=[];
    
    PUV_work.time = reshape(PUV_work.time, hourlen, segtotal);
    PUV_work.P = reshape(PUV_work.P, hourlen, segtotal);
    PUV_work.BuoyCoord.U = reshape(PUV_work.BuoyCoord.U, hourlen, segtotal);
    PUV_work.BuoyCoord.V = reshape(PUV_work.BuoyCoord.V, hourlen, segtotal);
    PUV_work.BuoyCoord.W = reshape(PUV_work.BuoyCoord.W, hourlen, segtotal);
    
    PUV_work.ShorenormalCoord.U = reshape(PUV_work.ShorenormalCoord.U, hourlen, segtotal);
    PUV_work.ShorenormalCoord.V = reshape(PUV_work.ShorenormalCoord.V, hourlen, segtotal);
    %% Detide PUV work
    disp('detide')
    for gg = 1:size(PUV_work.time,2)
        PUV_work.depth(gg) = nanmean(PUV_work.P(:,gg)) + doffp;
        x = PUV_work.P(:,gg);
        if isempty(find(isnan(x)))
            timedat = (1:1:size(PUV_work.time,1))'/2-1/2;
            tide = [ones(length(timedat),1), ...
            sin(2*pi*timedat/12.42/3600),cos(2*pi*timedat/12.42/3600), ... % M2
            sin(2*pi*timedat/12/3600),cos(2*pi*timedat/12/3600), ... % S2
            sin(2*pi*timedat/23.93/3600),cos(2*pi*timedat/23.93/3600)];% K1
            trend = tide * inv(tide.' * tide) * tide.' * x;
            PUV_work.P(:,gg) = x - trend;
        end

        x = PUV_work.ShorenormalCoord.U(:,gg);
        if isempty(find(isnan(x)))
            timedat = (1:1:size(PUV_work.time,1))'/2-1/2;
            tide = [ones(length(timedat),1), ...
            sin(2*pi*timedat/12.42/3600),cos(2*pi*timedat/12.42/3600), ... % M2
            sin(2*pi*timedat/12/3600),cos(2*pi*timedat/12/3600), ... % S2
            sin(2*pi*timedat/23.93/3600),cos(2*pi*timedat/23.93/3600)];% K1
            trend = tide * inv(tide.' * tide) * tide.' * x;
            PUV_work.ShorenormalCoord.U(:,gg) = x - trend;
        end

        x = PUV_work.ShorenormalCoord.V(:,gg);
        if isempty(find(isnan(x)))
            timedat = (1:1:size(PUV_work.time,1))'/2-1/2;
            tide = [ones(length(timedat),1), ...
            sin(2*pi*timedat/12.42/3600),cos(2*pi*timedat/12.42/3600), ... % M2
            sin(2*pi*timedat/12/3600),cos(2*pi*timedat/12/3600), ... % S2
            sin(2*pi*timedat/23.93/3600),cos(2*pi*timedat/23.93/3600)];% K1
            trend = tide * inv(tide.' * tide) * tide.' * x;
            PUV_work.ShorenormalCoord.V(:,gg) = x - trend;
        end
    end
    disp('Saving PUV_work')
    save([data_dir '/Level2_QC/' filename '_PUV_work.mat'],'PUV_work', '-v7.3')
    toc
    %% Compute different wave stats
    disp('Getting wave stats, writing to PUV_process')
    
    P = PUV_work.P;
    U = PUV_work.BuoyCoord.U;
    V = PUV_work.BuoyCoord.V;
    W = PUV_work.BuoyCoord.W;
    fs = PUV_work.fs;
    time = PUV_work.time;
    
    Uprime = PUV_work.ShorenormalCoord.U;
    Vprime = PUV_work.ShorenormalCoord.V;
    for gg = 1:size(P,2)
         if isempty(find(isnan(P(:,gg))))
            %% Wave statistics
            PUV_work_new.P(:,gg)=P(:,gg);
            PUV_work_new.time(:,gg)=time(:,gg);
            %[PUV_process(gg)] = vector_wave_stats_mspec_cutoff(Uprime(:,gg), Vprime(:,gg), P(:,gg),time(:,gg), doffp, doffu);
            [PUV_process(gg)] = vector_wave_stats_hanning(Uprime(:,gg), Vprime(:,gg), P(:,gg),time(:,gg), doffp, doffu, PUV_work.depth(gg));
         else
             PUV_process(gg).time = time(:,gg);
             PUV_process(gg).Spec = [];
             PUV_process(gg).ztest = [];
             PUV_process(gg).Hsig = [];
             PUV_process(gg).FC = [];
             PUV_process(gg).coh = [];
             PUV_process(gg).dir = [];
             PUV_process(gg).RS = [];
             PUV_process(gg).Eflux = [];
             PUV_process(gg).ids = [];
         end
    end
    %% Save data
    disp('saving')
    save([data_dir '/Level2_QC/' filename '_PUV_process_hanning.mat'], 'PUV_process', '-v7.3')
    toc

end

%% QC Level 2.5 - (NaN, Z-test, Reflection)
for ll = 1:size(deploy_notes,1)
    clear mm Ztest_ss Ztest_ig id R2_ss R2_ig
    filename = [char(string(deploy_notes.Location(ll))) '_' char(string(deploy_notes.Date(ll))) '_MOP' num2str(deploy_notes.MOP(ll)) '_' num2str(deploy_notes.Depth(ll)) 'm']
    load(fullfile(data_dir, 'Level2_QC', [filename '_PUV_process_hanning.mat']), 'PUV_process')
    load(fullfile(data_dir, 'Level2_QC', [filename '_PUV_work.mat']), 'PUV_work')
    
    disp('Remove PUV_process empty files')
    for mm = length(PUV_process):-1:1
        if isempty(PUV_process(mm).Spec)
            mm
            PUV_work.time(:,mm)=[];
            PUV_work.P(:,mm)=[];
            PUV_work.BuoyCoord.U(:,mm)=[];
            PUV_work.BuoyCoord.V(:,mm)=[];
            PUV_work.BuoyCoord.W(:,mm)=[];
            PUV_work.ShorenormalCoord.U(:,mm)=[];
            PUV_work.ShorenormalCoord.V(:,mm)=[];

            PUV_process(mm)=[];
        end
    end

    %%% Do Z-test check
    disp('Z-test')
    for mm = 1:length(PUV_process)
        Ztest_ss(mm) = PUV_process(mm).ztest.ztest_ss_sum;
        Ztest_ig(mm) = PUV_process(mm).ztest.ztest_ig_sum;
    end
    
    id = unique([find(Ztest_ig < 0.8) find(Ztest_ig > 1.2) find(Ztest_ss < 0.95) find(Ztest_ss > 1.05)]);
    PUV_work.time(:,id)=[];
    PUV_work.P(:,id)=[];
    PUV_work.BuoyCoord.U(:,id)=[];
    PUV_work.BuoyCoord.V(:,id)=[];
    PUV_work.BuoyCoord.W(:,id)=[];
    PUV_work.ShorenormalCoord.U(:,id)=[];
    PUV_work.ShorenormalCoord.V(:,id)=[];
    PUV_process(id) = []; 
    clear id

    %%% Reflection coefficient check
    disp('Reflection Coefficient Check')
    i_ig_sn = PUV_process(1).ids.i_ig;
    i_swell_sn = PUV_process(1).ids.i_swell;
    fm_sn = PUV_process(1).Spec.fm;
    df_sn = fm_sn(2)-fm_sn(1);

    for mm = 1:length(PUV_process)
        R2_ig(mm) = sum(PUV_process(mm).Eflux.negX(i_ig_sn).*df_sn)/sum(PUV_process(mm).Eflux.posX(i_ig_sn).*df_sn);
        R2_ss(mm) = sum(PUV_process(mm).Eflux.negX(i_swell_sn).*df_sn)/sum(PUV_process(mm).Eflux.posX(i_swell_sn).*df_sn);
        PUV_process(mm).MOP = deploy_notes.MOP(ll);
        PUV_process(mm).Location = deploy_notes.Location(ll);
        PUV_process(mm).Year = deploy_notes.Date(ll);
    end
    id = unique([find(R2_ss > 0.25) find(R2_ig > 2.5)]);
    PUV_work.time(:,id)=[];
    PUV_work.P(:,id)=[];
    PUV_work.BuoyCoord.U(:,id)=[];
    PUV_work.BuoyCoord.V(:,id)=[];
    PUV_work.BuoyCoord.W(:,id)=[];
    PUV_work.ShorenormalCoord.U(:,id)=[];
    PUV_work.ShorenormalCoord.V(:,id)=[];
    PUV_process(id) = []; 
    clear id

    %%% Make sure to include location and MOP

    PUV_process = orderfields(PUV_process, {'MOP', 'Location', 'Year','time', 'ids', 'Spec', 'Hsig', 'Eflux', 'FC', 'RS', 'coh', 'dir', 'ztest'});
    
    save(fullfile(data_dir, 'Level2_QC', 'temp', [filename '_PUV_process_hanning.mat']), 'PUV_process')
    save(fullfile(data_dir, 'Level2_QC', 'temp', [filename '_PUV_work_hanning.mat']), 'PUV_work')
    toc
end

%% Combine PUV_process, time, P, U and V
for ll = 1:size(deploy_notes,1)
    tic
    clearvars -except files deploy_notes names ll data_dir code_dir *_all time
    disp('New run')
    filename = [char(string(deploy_notes.Location(ll))) '_' char(string(deploy_notes.Date(ll))) '_MOP' num2str(deploy_notes.MOP(ll)) '_' num2str(deploy_notes.Depth(ll)) 'm']
    
    load(fullfile(data_dir, 'Level2_QC', 'temp', [filename '_PUV_work_hanning.mat']), 'PUV_work')
    load(fullfile(data_dir, 'Level2_QC', 'temp', [filename '_PUV_process_hanning.mat']), 'PUV_process')
    
    if ll == 1
        P_all = PUV_work.P;
        time = PUV_work.time;
        U_all = PUV_work.ShorenormalCoord.U;
        V_all = PUV_work.ShorenormalCoord.V;
        PUV_process_all = PUV_process;
    else
        P_all = [P_all PUV_work.P];
        time = [time PUV_work.time];
        U_all = [U_all PUV_work.ShorenormalCoord.U];
        V_all = [V_all PUV_work.ShorenormalCoord.V];
        PUV_process_all = [PUV_process_all PUV_process];
    end
end

id=[];
for ii = size(P_all,2):-1:1
    if ~isempty(find(isnan(P_all(:,ii))))
        id = [id ii];
    end

end

P_all(:,id)=[]; U_all(:,id)=[]; V_all(:,id)=[]; time(:, id)=[]; PUV_process_all(id)=[];
P = P_all; U = U_all; V = V_all; save(fullfile(data_dir, 'Level3_QC', 'PUV_work_P_hanning.mat'), 'P', 'U', 'V', 'time')
PUV_process = PUV_process_all; save(fullfile(data_dir, 'Level3_QC', 'PUV_process_hanning.mat'), 'PUV_process', '-v7.3')

%% QC Level 3 - PUV_incoming_outgoing
clearvars -except data_dir
load(fullfile(data_dir, 'Level3_QC', 'PUV_work_P_hanning.mat'), 'P', 'U', 'time')
load(fullfile(data_dir, 'Level3_QC', 'PUV_process_hanning.mat'), 'PUV_process')

dt = 0.5; fcutoff = 0.25;
for gg = 1:size(P,2)
    gg
    PUV(gg).time = time(:,gg);
    PUV_process(gg).depth = PUV_process(gg).Eflux.depth;

    [etaAa, ~, ~, ~]= etasolver_detide(P(:,gg), time(:,gg), size(P,1), PUV_process(gg).depth);
    [sigIN, sigOUT] = get_directionality_detrend(etaAa,U(:,gg),'shoreward', PUV_process(gg).depth, dt);

    PUV(gg).eta_IN = sigIN; 
    [~, etaAss, etaAig]= timeseries_ss_ig(sigIN, [0:0.5:1]);

    PUV(gg).eta_IN_ss = etaAss;
    PUV(gg).eta_IN_ig = etaAig;

    eta = detrend(sigIN); if isrow(eta); eta = eta';end
    [fm, Spp, ~, ~, ~, ~] = get_spectrum(eta, 7200, 2, 0.05);
    ii = find(fm <= fcutoff); id_low = find(min(abs(fm - 0.004))== abs(fm-0.004));
    Spp(1:id_low) = NaN; SSEt = Spp(ii); fmt = fm(ii);
    PUV_process(gg).Spec.Seta_in = SSEt;

    eta = detrend(sigOUT); if isrow(eta); eta = eta';end
    [fm, Spp, ~, ~, ~, ~] = get_spectrum(eta, 7200, 2, 0.05);
    ii = find(fm <= fcutoff); id_low = find(min(abs(fm - 0.004))== abs(fm-0.004));
    Spp(1:id_low) = NaN; SSEt = Spp(ii); fmt = fm(ii);
    PUV_process(gg).Spec.Seta_out = SSEt;

    [t_verified,~,verified,~,~] = getNOAAtide(time(1,gg),time(end,gg),'9410230');
    PUV_process(gg).tide.t = t_verified;
    PUV_process(gg).tide.z = verified;
    PUV_process(gg).tide.depth = PUV_process(gg).depth - mean(verified);
end
save(fullfile(data_dir, 'Level3_QC', 'PUV_process_hanning.mat'), 'PUV_process', '-v7.3')
save(fullfile(data_dir, 'Level3_QC', 'PUV_eta.mat'), 'PUV', '-v7.3')


%% QC Level 3 - Moments, Boundwave & bispectra
load(fullfile(data_dir, 'Level3_QC', 'PUV_process_hanning.mat'), 'PUV_process')
load(fullfile(data_dir, 'Level3_QC', 'PUV_eta.mat'), 'PUV')
fs = 2;
for gg = 1:length(PUV)
    gg
    % Moments
    disp('Computing moments')
    inorm = 1;
    T = length(PUV(gg).eta_IN)/fs;
    flim = [0.004 0.04 0.25];
    [m2_in,m3_in] = sectormoments(T, PUV(gg).eta_IN, flim, inorm);
    moments(gg).m2_in = m2_in;
    moments(gg).m3_in = m3_in;

    % Boundwave
    disp('Computing 1D boundwave')
    % 1D forced boundwave
    [Z_bound_in] = boundwave_zig1D(PUV(gg).eta_IN, PUV_process(gg).depth);
    boundwave(gg).Z_1Dbound_in = Z_bound_in;

    disp('Computing 2D boundwave')
    % 2D boundwave
    [Z_IG, E_bound, f_bound, alpha] = boundwave_zig2D(PUV(gg).eta_IN, PUV_process(gg).depth, PUV_process(gg).FC, PUV_process(gg).Spec.fm, 'PUV');
    boundwave(gg).Z_2Dbound_in = Z_IG;
    boundwave(gg).E_bound = E_bound;
    boundwave(gg).f_bound = f_bound;
    boundwave(gg).alpha = alpha;

    % Bispectra
    disp('Computing bispectra')
    [data] = bispec(PUV(gg).eta_IN, boundwave(gg).alpha, boundwave(gg).f_bound);
    bispec_var(gg).var_alpha = data;
    [data] = bispec(PUV(gg).eta_IN, ones(size(boundwave(gg).alpha,1), size(boundwave(gg).alpha,2)), boundwave(gg).f_bound);
    bispec_var(gg).const_alpha = data;

end
% Save data

disp('saving')
save(fullfile(data_dir, 'Level3_QC', 'PUV_moments_hanning.mat'), 'moments', '-v7.3')
save(fullfile(data_dir, 'Level3_QC', 'PUV_boundwave_hanning.mat'), 'boundwave', '-v7.3')
save(fullfile(data_dir, 'Level3_QC', 'temp', 'PUV_bispec_hanning.mat'), 'bispec_var', '-v7.3')
        


%%
for gg = 1001%:length(PUV_process)
    clearvars -except PUV PUV_process boundwave gg data_dir bispec_var
    gg
    tic
    disp('Computing 2D boundwave')
    % 2D boundwave
    [Z_IG, E_bound, f_bound, alpha] = boundwave_zig2D(PUV(gg).eta_IN, PUV_process(gg).depth, PUV_process(gg).FC, PUV_process(gg).Spec.fm, 'PUV');
    boundwave(gg).Z_2Dbound_in = Z_IG;
    boundwave(gg).E_bound = E_bound;
    boundwave(gg).f_bound = f_bound;
    boundwave(gg).alpha = alpha;
    toc
    disp('Computing bispectra')
    [data] = bispec(PUV(gg).eta_IN, boundwave(gg).alpha, boundwave(gg).f_bound);
    bispec_var(gg).var_alpha = data;
    [data] = bispec(PUV(gg).eta_IN, ones(size(boundwave(gg).alpha,1), size(boundwave(gg).alpha,2)), boundwave(gg).f_bound);
    bispec_var(gg).const_alpha = data;
    toc
    if rem(gg,50) == 0
        save(fullfile(data_dir, 'Level3_QC', 'temp', ['PUV_bispec_' char(string(gg)) '_hanning.mat']), 'bispec_var', '-v7.3')
        save(fullfile(data_dir, 'Level3_QC', 'temp', ['PUV_boundwave_', char(string(gg)), '_hanning_incoming.mat']), 'boundwave', '-v7.3')
    end

end
% Save data
toc
disp('saving')
%save(fullfile(data_dir, 'Level3_QC', 'PUV_boundwave_hanning.mat'), 'boundwave', '-v7.3')
%save(fullfile(data_dir, 'Level3_QC', ['PUV_bispec_' char(string(gg)) '_hanning.mat']), 'bispec_var', '-v7.3')
toc


%% Pull Bulk statistics

for ii = 1:length(PUV_process)
    Bulk.time(ii) = PUV_process(ii).time(1);
    Bulk.depth(ii) = PUV_process(ii).depth;
    Bulk.MOP(ii)=PUV_process(ii).MOP;
    Bulk.tide(ii) = PUV_process(ii).tide;
    fx = gradient(smoothdata(Bulk.tide(ii).z, 'Gaussian', 10), Bulk.tide(ii).t);

    if all(fx < 0) % all neg - ebb
        Bulk.tide_lvl_description{ii} = 'ebb';
        Bulk.tide_lvl(ii) = 1;
    elseif all(fx > 0) % all pos - flood
        Bulk.tide_lvl_description{ii} = 'flood';
        Bulk.tide_lvl(ii) = 3;
    elseif fx(1) > 0 & fx(end) < 0 % high tide
        Bulk.tide_lvl_description{ii} = 'high';
        Bulk.tide_lvl(ii) = 4;
    elseif fx(1) < 0 & fx(end) > 0 % low tide
        Bulk.tide_lvl_description{ii} = 'low';
        Bulk.tide_lvl(ii) = 2;
    else
        Bulk.tide_lvl_description{ii} = 'unknown';
        Bulk.tide_lvl(ii) = 0;
    end

    i_ig_sn = PUV_process(ii).ids.i_ig;
    i_swell_sn = PUV_process(ii).ids.i_swell;
    fm_sn = PUV_process(ii).Spec.fm;
    df_sn = fm_sn(2)-fm_sn(1);

    Bulk.df = df_sn;
  
    Bulk.Hs_ss(ii) = 4 * sqrt( sum( PUV_process(ii).Spec.Seta_in(i_swell_sn) * Bulk.df ) );
    Bulk.Hs_ig(ii) = 4 * sqrt( sum( PUV_process(ii).Spec.Seta_in(i_ig_sn) * Bulk.df ) );
    Bulk.Hs_ig_out(ii) = 4 * sqrt( sum( PUV_process(ii).Spec.Seta_out(i_ig_sn) * Bulk.df ) );
    Bulk.Hs(ii) = 4 * sqrt( nansum( PUV_process(ii).Spec.Seta_in * Bulk.df ) );

    [fc, fsp, fkurt, fp] = get_fstats(PUV_process(ii).Spec.fm, PUV_process(ii).Spec.Seta_in, Bulk.df);
    Bulk.fc_in(ii)=fc;
    Bulk.fspread_in(ii)=fsp;
    Bulk.fpeak_in(ii)=fp;
    Bulk.Tm(ii) = PUV_process(ii).ids.Tm;
    Bulk.fpeak(ii) = PUV_process(ii).ids.fpeak;
    Bulk.fspread(ii) = PUV_process(ii).ids.fspread;
    Bulk.dir1_ss_sum(ii) = PUV_process(ii).dir.dir1_ss_sum;
    Bulk.dir2_ss_sum(ii) = PUV_process(ii).dir.dir2_ss_sum;
    Bulk.spread1_ss_sum(ii) = PUV_process(ii).dir.spread1_ss_sum;
    Bulk.spread2_ss_sum(ii) = PUV_process(ii).dir.spread2_ss_sum;
    Bulk.spread1_ss(ii) = PUV_process(ii).dir.spread1_ss;
    Bulk.spread2_ss(ii) = PUV_process(ii).dir.spread2_ss;
    Bulk.spread1_ig(ii) = PUV_process(ii).dir.spread1_ig;
    Bulk.spread2_ig(ii) = PUV_process(ii).dir.spread2_ig;

    a = Bulk.Hs_ss(ii)/2;
    d = PUV_process(ii).Eflux.depth;
    [k,omega]=dispersion(d, PUV_process(ii).Spec.fm');
    k=k(find(max(PUV_process(ii).Spec.Seta_in)==PUV_process(ii).Spec.Seta_in));
    Bulk.Ur(ii)=(a/d)./((k*d)^2);

    Bulk.R2_ig_sn(ii) = sum(PUV_process(ii).Eflux.negX(i_ig_sn).*df_sn) / sum(PUV_process(ii).Eflux.posX(i_ig_sn).*df_sn);
    Bulk.R2_ss_sn(ii) = sum(PUV_process(ii).Eflux.negX(i_swell_sn).*df_sn) / sum(PUV_process(ii).Eflux.posX(i_swell_sn).*df_sn);
    Bulk.mean_tide(ii) = nanmean(Bulk.tide(ii).z);
  
    Bulk_moments.bp_ig(ii) = rad2deg(atan2(imag(moments(ii).m3.dif),real(moments(ii).m3.dif)));
    Bulk_moments.a_ig(ii) = imag(moments(ii).m3.dif);
    Bulk_moments.s_ig(ii) = real(moments(ii).m3.dif);
    Bulk_moments.bc_ig(ii) = sqrt(dot(moments(ii).m3.dif,moments(ii).m3.dif,3)); if Bulk_moments.bc_ig(ii) < 0; Bulk_moments.bc_ig(ii) = Bulk_moments.bc_ig(ii) + 360;end
    Bulk_moments.bp_ss(ii) = rad2deg(atan2(imag(moments(ii).m3.self(2)),real(moments(ii).m3.self(2))));
    Bulk_moments.a_ss(ii) = imag(moments(ii).m3.self(2));
    Bulk_moments.s_ss(ii) = real(moments(ii).m3.self(2));
    Bulk_moments.bc_ss(ii) = sqrt(dot(moments(ii).m3.self(2), moments(ii).m3.self(2), 3));
end

for gg = 1:length(PUV_process)

    IG(gg).MOP = PUV_process(gg).MOP;
    IG(gg).Location = PUV_process(gg).Location;
    IG(gg).Year = PUV_process(gg).Year;
    IG(gg).time = PUV_process(gg).time(1);
    IG(gg).f = PUV_process(gg).Spec.fm;
    i_swell = find(IG(gg).f < 0.25 & IG(gg).f > 0.04);
    i_ig = find(IG(gg).f < 0.04 & IG(gg).f > 0.004);
    IG(gg).f_ss = PUV_process(gg).Spec.fm(i_swell);
    IG(gg).f_ig = PUV_process(gg).Spec.fm(i_ig);
    IG(gg).obs_in = smoothdata(PUV_process(gg).Spec.Seta_in, 'Gaussian', 25);
    IG(gg).obs_in_ss = IG(gg).obs_in(i_swell);
    IG(gg).obs_in_ig = IG(gg).obs_in(i_ig);

    IG(gg).obs_out = smoothdata(PUV_process(gg).Spec.Seta_out, 'Gaussian', 25);
    IG(gg).obs_out_ss = IG(gg).obs_out(i_swell);
    IG(gg).obs_out_ig = IG(gg).obs_out(i_ig);
    
    IG(gg).bound_ig = interp1(boundwave(gg).f_bound, boundwave(gg).E_bound, IG(gg).f_ig);
    IG(gg).free_ig = IG(gg).obs_ig - IG(gg).bound_ig;
    for ii = 1:length(IG(gg).free_ig)
        if IG(gg).free_ig(ii)< 0
            IG(gg).free_ig(ii)=NaN;
        end
    end
    
    IG(gg).frac_bound = sum(IG(gg).bound_ig)./sum(IG(gg).obs_ig);
    IG(gg).frac_free = sum(IG(gg).free_ig)./sum(IG(gg).obs_ig);
    frac_bound(gg) = sum(IG(gg).bound_ig)./sum(IG(gg).obs_ig);
    frac_free(gg) = sum(IG(gg).free_ig)./sum(IG(gg).obs_ig);
end

% Bulk.frac_bound = frac_bound;
% Bulk.frac_free = frac_free;

save(fullfile(data_dir, 'Level3_QC', 'Wave_stats_hanning_incoming.mat'), 'Bulk', 'Bulk_moments', 'IG')

