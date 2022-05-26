clear all
tic
%% Load all instrument data - need to create excel sheet from Reefbreak deployment sheet
% input directory where data is - one deloyment per folder
cd /Volumes/LANGE_Passport/PUV/Ruby2D/
files = dir(pwd)
for ii = 1:length(files);names{ii}=files(ii).name;end
deploy_notes = readtable(string(names(contains(names, 'Notes'))))

[~,ii]=sort(deploy_notes.MOP);
deploy_notes=deploy_notes(ii,:);
aa=tabulate(deploy_notes.MOP); 

if length(find(aa(:,2)>1))>0
    mm = find(aa(:,2)>1); % if multiple sensors on single mop line
    for jj = 1:length(mm)
        ll = mm(jj);
        nn = find(aa(ll,1)==deploy_notes.MOP);
        [~,ii]=sort(deploy_notes.Depth(nn));
        deploy_notes(nn,:)=deploy_notes(nn(ii),:);
    end
end
% remove any non-folder files
for nn = length(files):-1:1 
    if files(nn).name(1) == '.' | files(nn).isdir == 0
        files(nn)=[];
    end
end

%% Raw Process Data - Level 1
for ii = 1:length(files)
    directory = ['MOP' num2str(deploy_notes.MOP(ii)) '-' num2str(deploy_notes.Depth(ii)) 'm'];
    filename = ['Torrey_Ruby2D_' num2str(deploy_notes.MOP(ii)) '_' num2str(deploy_notes.Depth(ii)) 'm'];
    % input lat/lon of sensor
    LATLON = [deploy_notes.Latitude(ii) deploy_notes.Longitude(ii)];
    % input compass heading of beam 1 of ADV
    rot_angle = [deploy_notes.Heading(ii)];
    % input clock drift between beginning and end of deployment
    clockdrift = [deploy_notes.ClockDrift(ii)];
    [PUV] = PUV_raw_process(directory, filename, LATLON, rot_angle, clockdrift)

end
%% Compute Stats - Level 2

for ii = 1%:length(files)
    directory = ['MOP' num2str(deploy_notes.MOP(ii)) '-' num2str(deploy_notes.Depth(ii)) 'm'];
    filename = ['Torrey_Ruby2D_' num2str(deploy_notes.MOP(ii)) '_' num2str(deploy_notes.Depth(ii)) 'm'];
    load(['Level1_QC/' filename '_processed'], 'PUV')
    PUV.mop = filename(15:17);
    doffp = deploy_notes.Doffp(ii)/100;

    %% Rotate to shorenormal coordinates (MOP)
    
    [PUV.ShorenormalCoord.U, PUV.ShorenormalCoord.V, shorenormal] = rotate_shorenormal(PUV.BuoyCoord.U, PUV.BuoyCoord.V, PUV.mop);
    
    tic
    % 66cm sand to top of pressure port on deployment. 78cm on recovery. 
    % doffp = 0.76; % splitting the difference between the two
    % take doffp from Deployment Notes
    doffu = doffp+0.371; %adding 37.1 cm to u,v sample point
    PUV_work = PUV;
    
    %% Interp NaNs if less than 2sec gap
    PUV_work.P = fillmissing(PUV.P,'linear','SamplePoints',PUV.time,'MaxGap',seconds(2));
    PUV_work.BuoyCoord.U = fillmissing(PUV.BuoyCoord.U,'linear','SamplePoints',PUV.time,'MaxGap',seconds(2));
    PUV_work.BuoyCoord.V = fillmissing(PUV.BuoyCoord.V,'linear','SamplePoints',PUV.time,'MaxGap',seconds(2));
    PUV_work.BuoyCoord.W = fillmissing(PUV.BuoyCoord.W,'linear','SamplePoints',PUV.time,'MaxGap',seconds(2));
    PUV_work.ShorenormalCoord.U = fillmissing(PUV.ShorenormalCoord.U,'linear','SamplePoints',PUV.time,'MaxGap',seconds(2));
    PUV_work.ShorenormalCoord.V = fillmissing(PUV.ShorenormalCoord.V,'linear','SamplePoints',PUV.time,'MaxGap',seconds(2));
    PUV_work.T = fillmissing(PUV.T,'linear','SamplePoints',PUV.time,'MaxGap',seconds(2));
    toc
    %% Split data into hour segments
    
    for ii = 0:floor(length(PUV_work.time)/7200)
        if minute(PUV.time(1+7200*ii)) ~= 0
            PUV_work.time(1+7200*ii)=[];
            PUV_work.P(1+7200*ii)=[];
            PUV_work.BuoyCoord.U(1+7200*ii)=[];
            PUV_work.BuoyCoord.V(1+7200*ii)=[];
            PUV_work.BuoyCoord.W(1+7200*ii)=[];
            PUV_work.ShorenormalCoord.U(1+7200*ii)=[];
            PUV_work.ShorenormalCoord.V(1+7200*ii)=[];
        end
    end

    hourlen = (60*60*PUV_work.fs);
    segtotal = floor(length(PUV_work.time)/hourlen);
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

    
    PUV = PUV_work;
    save(['Level2_QC/' filename '_PUV_work'],'PUV', '-v7.3')
    toc
    %% Compute different wave stats
    disp('Getting wave stats, writing to PUV_process')
    tic
    P = PUV_work.P;
    U = PUV_work.BuoyCoord.U;
    V = PUV_work.BuoyCoord.V;
    W = PUV_work.BuoyCoord.W;
    fs = PUV_work.fs;
    time = PUV_work.time;
    toc
    
    %% Shorenormal, +x onshore, +y alongshore north, +z upward 
    % use for boundwave codes
    Uprime = PUV_work.ShorenormalCoord.U;
    Vprime = PUV_work.ShorenormalCoord.V;
    for ii = 1:size(P,2)
        if isempty(find(isnan(P(:,ii))))
%             [etaAa, etaAss, etaAig, tm]= etasolver(P(:,ii), time(:,ii), size(P,1));
%             eta_level(ii).P = P(:,ii);
%             eta_level(ii).etaAa = etaAa;
%             eta_level(ii).etaAss = etaAss;
%             eta_level(ii).etaAig = etaAig;
%             eta_level(ii).time = time;
            %[eta_level(ii).m2,eta_level(ii).m3] = sectormoments(7200/2, eta_level(ii).P , [0.004 0.04 0.25],1);
            [PUV_process(ii)] = vector_wave_stats_mspec(Uprime(:,ii), Vprime(:,ii), P(:,ii),time(:,ii), doffp, doffu);
        
        end
    end
    
    % Buoy coordinates: +x West, +y North, +z Up
    Uprime = PUV_work.BuoyCoord.U;
    Vprime = PUV_work.BuoyCoord.V;
    for ii = 1:size(P,2)
        if isempty(find(isnan(P(:,ii))))
%             [etaAa, etaAss, etaAig, tm]= etasolver(P(:,ii), time(:,ii), size(P,1));
%             eta_level(ii).P = P(:,ii);
%             eta_level(ii).etaAa = etaAa;
%             eta_level(ii).etaAss = etaAss;
%             eta_level(ii).etaAig = etaAig;
%             eta_level(ii).time = time;
            %[eta_level(ii).m2,eta_level(ii).m3] = sectormoments(7200/2, eta_level(ii).P , [0.004 0.04 0.25],1);
            [PUV_process_buoy(ii)] = vector_wave_stats_mspec(Uprime(:,ii), Vprime(:,ii), P(:,ii),time(:,ii), doffp, doffu);
        
        end
    end
    toc

    %save('eta', 'eta_level', '-v7.3')
    save(['Level2_QC/' filename '_PUV_process.mat'], 'PUV_process', 'PUV_process_buoy', '-v7.3')

toc

end
