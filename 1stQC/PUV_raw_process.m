% PRE-PROCESS VECTOR DATA
% Load Nortek Vector files
% Account for timedrift, battery failure, beginning of deployment & inspection
% removing bad data (minCorr < 70%)
% rotate to buoy coordinates (+x WEST, +y NORTH, +z up)
% save to PUV variable
% PUV =     
%       time (datetime)
%       fs (Sampling Frequency Hz)
%       P (Pressure dBar)
%       U (rotated U +x WEST)
%       V (rotate V +y NORTH)
%       W (+z up)
%       Vel.U & Vel.V (XYZ coordinate velocities : +x crosshore, onshore, +y alongshore, south)
%       LATLON ([LAT LON])
%       rotation.sensor (clockwise angle from magnetic north or direction of probe 1, X)
%       rotation.mag (clockwise angle of magnetic north from true north, calculated at the beginning of deployment)
%       inttime.start / inttime.end (start & end times of data segments, NaN's inbetween, datetime)
%

%%
function [PUV] = PUV_raw_process(directory, filename, LATLON, rot_angle, clockdrift)%,inspec_id)
    
    %% ----------------- GEOGRAPHICAL INFO -----------------

    LAT = LATLON(1);  LON = LATLON(2);

    %% ----------------- File location -----------------

    % Determine beginning string name from files
    filenames = struct2cell(dir([directory '/*.dat'])); filenames = filenames(1,:);
    depstr = filenames{1}(1:find(any(diff(char(filenames(:))),1),1,'first')-1)

    % load SEN and DAT file
    % input beginning and end dates
    for ii = 1:length(filenames)
        id = char(string(ii));
        senfile = sprintf('%s/%s%s.sen',directory,depstr,id)
        datfile = sprintf('%s/%s%s.dat',directory,depstr,id)

        eval(['SEN' id ' = load(senfile); DAT' id ' = load(datfile);'])

        eval(['date(ii,:) = [SEN' id '(1,3) SEN' id '(1,1:2) SEN' id '(1,4:6)]'])
        eval(['date_end(ii,:) = [SEN' id '(end,3) SEN' id '(end,1:2) SEN' id '(end,4:6)]'])
    end
    
    
    
    
    %% ----------------- Time drift -----------------
    % find the sampling rate from hdr file (NEED TO ADD 1Hz option in step above)
    fid=fopen(sprintf('%s/%s%s.hdr',directory,depstr,'1')); 
    linenum = 12;
    C = split(string(textscan(fid,'%s',1,'delimiter','\n', 'headerlines',linenum-1)));

    if isempty(find(C=='Sampling'))
        sprintf('ERROR: Wrong line for Sampling Rate')
    else
        fs = double(C(3)) % save sampling rate
    end
    
    % find total duration of deployment (including intermitent bursts due to battery depletion 
    totalseconds = seconds(time(caldiff([datetime(date(1,:), 'InputFormat', 'YYYYMMDDHHmmSS') datetime(date_end(end,:), 'InputFormat', 'YYYYMMDDHHmmSS')+ seconds(1/fs)], 'Time')))
    % create linspace that is length of deployment, and starts at 0 to the drifttime
    driftfix = linspace(seconds(0), seconds(clockdrift), totalseconds*fs+1);
   %% ----------------- Save original data to netcdf -----------------
    save_raw_data_netcdf
    %% ----------------- Find beginning of intermintent data recording -----------------
    % Assuming a 2Hz sampling scheme - SEN file should have a value
    % recorded every second. If the difference between neighbooring values >
    % 1s, then that is where the data is truncated to.

    cutoff = [];
    for ii = 1:length(filenames)
        id=char(string(ii))
        % generate datetime array from SEN recorded values
        eval(['date_temp= datetime([SEN' id '(:,3) SEN' id '(:,1:2) SEN' id '(:,4:6)], ''InputFormat'', ''YYYYMMDDHHmmSS'');']);
        
        % check if timestamps are further than 1s apart -> if so, that's
        % the beginning of cutoff
        aa = find(diff(date_temp) > seconds(1));
        if ~isempty(aa)
            cutoff.burst = ii; 
            cutoff.id = aa(1);
        end
        clear date_temp
        
        if ~isempty(cutoff)
            sprintf('Stopping timeseries at %s', id)
            return
        end
    end

    
    % If a cutoff was found: Looking in that set of data, remove all values
    % past the cutoff time, and fix the date_end to be the new end time.
    % Also delete any files that might exist past that cutoff. 

    if ~isempty(cutoff)
        id = char(string(cutoff.burst))
        
        % delete any values past cutoff in current interval
        % need to make DAT file as long as double the SEN file (assuming 2Hz sampling)
        eval(['DAT' id '(end+1:end+(length(SEN' id ')*fs-length(DAT' id ')),:) = NaN;'])
        eval(['SEN' id '(cutoff.id:end,:)=[];'])
        eval(['DAT' id '(cutoff.id*fs-1:end,:)=[];'])
        eval(['date_end(' id ',:) = [SEN' id '(end,3) SEN' id '(end,1:2) SEN' id '(end,4:6)]'])

        % delete any data intervals past cutoff
        for ii = cutoff.burst+1:length(filenames)   
            id = char(string(ii))
            eval(['clear SEN' id ' DAT' id '; date(' id ',:) = []; date_end(' id ',:)=[];'])
        end
    end

    %% ----------------- Combine all the data ----------------- 
    
    % make date & date_end into datetime format
    date = datetime(date, 'InputFormat', 'YYYYMMDDHHmmSS')
    date_end = datetime(date_end, 'InputFormat', 'YYYYMMDDHHmmSS')
     
    % create datetime array from start to end of deployment
    full_date = date(1):seconds(1/fs):date_end(end)+seconds(1/fs);
        
    DAT = NaN(length(full_date), size(DAT1,2));
    SEN = NaN(length(full_date), size(SEN1,2));
    
    if fs == 2
        sprintf('Sampling: 2Hz')
        for ii = 1:cutoff.burst
            id=char(string(ii))
            % add NaN's to end of DAT file to match 2*length of SEN file. 
            eval(['DAT' id '(end+1:end+(length(SEN' id ')*2-length(DAT' id ')),:) = NaN;'])
            % add DAT data to new array for full deployment
            eval(['DAT(find(date(ii) == full_date):find(date_end(ii) == full_date)+1,:) = DAT' id ';'])
            
            eval(['SEN(find(date(ii) == full_date):2:find(date_end(ii) == full_date)+1,:) = SEN' id ';'])
            eval(['SEN(find(date(ii) == full_date)+1:2:find(date_end(ii) == full_date),:) = SEN' id '(1:end-1,:);'])
        end
    else fs == 1
        sprintf('Sampling: 1Hz')
        for ii = 1:cutoff.burst
            id=char(string(ii))
            % add NaN's to end of DAT file to match 2*length of SEN file. 
            eval(['DAT' id '(end+1:end+(length(SEN' id ')-length(DAT' id ')),:) = NaN;'])
            % add DAT data to new array for full deployment
            eval(['DAT(find(date(ii) == full_date):find(date_end(ii) == full_date)+1,:) = DAT' id ';'])
            
            eval(['SEN(find(date(ii) == full_date):find(date_end(ii) == full_date)+1,:) = SEN' id ';'])
        end
    end

    % adjusting the drift array to new length of data
    driftfix(length(full_date)+1:end)=[];
    % adjusting for drift
    full_date = full_date - driftfix; 
    
    %% TODO AUTOMATE ----------------- Remove above surface level points & Heading/Tilt issues -----------------
    % need to fix to get input
    
    % There shouldn't be any crazy changes in Heading, Pitch or Roll
    % (no greater than 10deg) (if between 10 and 20deg, should do 'map to
    % vertical', but here those instances agree with deploy and inspection
    % - issues with P, so we replace with NaN). 
    
    % P should have dips to 0 (at surface), at the beginning of the
    % deployment and during any inspections. Look at the excel sheet to
    % check where there should be confirmed dips, and replace with NaNs. 

%     figure(1);clf;
%     subplot(221)
%     plot(DAT(:,15));
%     title('P')
%     subplot(222)
%     plot(SEN(:,11));
%     title('Heading')
%     subplot(223)
%     plot(SEN(:,12));
%     title('Pitch')
%     subplot(224)
%     plot(SEN(:,13));
%     title('Roll')
%     
    depl_start_id = 223649;
    inspe_start_id = 13123100;%13345000;
    inspe_end_id = 13126300;%13345600;
    pitch_issue_start_id = 13343100;
    pitch_issue_end_id = 13346400;
    
    inspec_id = [1 depl_start_id; 
              inspe_start_id inspe_end_id;
              pitch_issue_start_id pitch_issue_end_id];

    for ii = size(inspec_id,1):-1:2
        full_date(inspec_id(ii,1):inspec_id(ii,2))=NaT;
        DAT(inspec_id(ii,1):inspec_id(ii,2),:)=NaN;
        SEN(inspec_id(ii,1):inspec_id(ii,2),:)=NaN; 
    end
    
    full_date(inspec_id(1,1):inspec_id(1,2))=[];
    DAT(inspec_id(1,1):inspec_id(1,2),:)=[];
    SEN(inspec_id(1,1):inspec_id(1,2),:)=[]; 
    
    figure(2);clf;
    subplot(221)
    plot(DAT(:,15));
    title('P')
    subplot(222)
    plot(SEN(:,11));
    title('Heading')
    subplot(223)
    plot(SEN(:,12));
    title('Pitch')
    subplot(224)
    plot(SEN(:,13));
    title('Roll')
    
    %% ----------------- Remove all values with minCorr < 70% -----------------
    bad_data = find(min(DAT(:,[12:14]),[],2) < 70);

    full_date(bad_data) = NaT;
    DAT(bad_data, :) = NaN;
    SEN(bad_data, :) = NaN;

    %% Pull data into separate variables

    U = DAT(:,3); V = DAT(:,4); W = DAT(:,5); % in m/s
    T = SEN(:, 14);
    P = DAT(:,15);

    %% Rotate to true North - can then rotate to shorenormal MOP
    % Compute the magnetic declination based on lat/lon and start time of
    % deployment when compass heading was taken). 

    % Check that ADV set to XYZ coordinate system. 
    % Assume ADV beam 1 pointing to shore, coordinate system is +x onshore, and
    % +y alongshore toward the south. 

    % Goal is to get velocity values into the coordinate system of MOPS (+x
    % WEST, +y NORTH). So +x is going to 270deg, but because +x in the
    % original POV is the compass heading + magentic declination off of NORTH,
    % subtract that from 270 to get clockwise rotation values (need to flip
    % sign in rotation matrix). 

    % if sensor upward looking: +z is down

    % ALL THIS SHOULD BE CONFIRMED BY COMPARING FC WITH BUOY


    [~,~, magDeclination, ~,~] = igrfmagm(0, LAT,LON, decyear(SEN1(1,3),SEN1(1,1),SEN1(1,2)), 13)

    fid=fopen(sprintf('%s/%s%s.hdr',directory,depstr,'1')); 
    linenum = 32;
    C = split(string(textscan(fid,'%s',1,'delimiter','\n', 'headerlines',linenum-1)));

    if isempty(find(C=='Coordinate'))
        sprintf('ERROR: Wrong line for Coordinate System')
    else
        sprintf('Coordinate System = %s', C(3))
    end

    if C(3)== 'XYZ'

        theta1_mag = rot_angle; % CHECK THAT IT SHOULD BE PULLED DIRECTLY FROM VALUE
        theta1_true = theta1_mag+ magDeclination; % XYZ coordinates at deg true

        rotation = 270 - theta1_true; % +x is WEST, +y is NORTH
        alpha_rad = rotation*pi/180; 
        rot_mat   = [cos(alpha_rad) sin(alpha_rad);
                 -sin(alpha_rad)  cos(alpha_rad)];

        uv=rot_mat*[U';V'];
        U_true=uv(1,:).';       % this should be in true North coordinates now! 
        V_true=uv(2,:).';
    elseif C(3) == 'ENU'
        printf('Does a coordinate system rotation need to happen?')
    end

   

    % NOW: +x is WEST, +y is NORTH
    %%
    
    %% SAVE VARIABLES

    PUV.time = full_date'; 
    PUV.fs = fs;
    PUV.P = P;
    PUV.BuoyCoord.U = U_true;
    PUV.BuoyCoord.V = V_true;
    PUV.BuoyCoord.W = -W;
    PUV.InstrCoord.U = U;
    PUV.InstrCoord.V = V;
    PUV.T = T;
    PUV.LATLON = LATLON;
    PUV.rotation.sensor = rot_angle;
    PUV.rotation.mag = magDeclination;
    

    hourlen = (60*60*PUV.fs)
    segtotal = floor(length(PUV.time)/hourlen)
    aa = find(hour(PUV.time(1:hourlen)) == hour(PUV.time((1))));
    if size(aa) < hourlen
        sprintf('Beginning timeseries at the next hour')
        PUV.P(aa)=[];
        PUV.BuoyCoord.U(aa)=[];
        PUV.BuoyCoord.V(aa)=[];
        PUV.BuoyCoord.W(aa)=[];
        PUV.InstrCoord.U(aa)=[];
        PUV.InstrCoord.V(aa)=[];
        PUV.T(aa)=[];
        PUV.time(aa)=[];

    end
    

    save([filename '_processed'], 'PUV')
    save_raw_data_netcdf
end
