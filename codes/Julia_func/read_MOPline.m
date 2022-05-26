function MOP = read_MOPline(MOPname,date1,date2)
% MOP = readMOPline(MOPname,date1,date2)
% This function reads the 1D E(f) spectrum from the CDIP thredds server
% plus some other fun stuff. It can be altered to retrieve other stuff as
% well, I just care about 1D spectra at the moment.
% Based on Bonnie Ludka's getwaves function.
%
% INPUT:
% MOPname: string format, example: 'D0667' for Cardiff 2012 experiment
% date1, date2: bounding dates in datetime format, ex: datetime(2010,1,1);
%
% OUTPUT:
% MOP structure containing
%   time:         in datetime format, hourly
%   spec1D:       a 1D E(f) spectra, dimensions: (time,frequency)
%   frequency:    vector 0.04<f<0.4
%   Hs:           significant wave height, dimensions: time
%   Dp:           peak wave direction, in degrees, dimensions: time  
%   depth:        depth of MOP line
%   lat:          latitude  
%   lon:          longitude  
%   shorenormal:  direction in degrees (0-360) of shore normal  
%   name:         MOP name, should be equal to MOPname input
%   flag1:        QC flag (primary)
%                 1:good, 2: not evaluated, 3: questionable,
%                 4: bad, 5: missing
%   flag2:        QC flag (secondary)
%                 0: unspecified, 1: insufficient input, 2: low energy
%   fbounds:      lower and upper bounds of frequency bins
%   a1,b1,a2,b2:  first 4 fourier components, used for directional
%                 information and MEM 2D spectrum estimations
%
% DERIVED QTYS:
%   fp:           1/peak period, in Hz, dimensions: time
%   dirspread1:   directional spread 1, in degrees (closest to line
%                 definition of standard deviation)
%   dirspread2:   directional spread 2, in degrees (use if reflection is
%                 important!)
%
% For more info, see http://cdip.ucsd.edu/documentation

% Julia Fiedler, jfiedler@ucsd.edu March 2019


% TODO: nowcast, hindcast!!
tSwitch = datetime(2018,1,1,0,0,0);

if date2<tSwitch
    url = ['http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/model/MOP_alongshore/' MOPname '_hindcast.nc'];
    disp('Using hindcast')
elseif date1>=tSwitch
    if date1 == dateshift(datetime,'start','day')
    url = ['http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/model/MOP_alongshore/' MOPname '_forecast.nc'];
    disp('Using forecast')
    else
    
    url = ['http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/model/MOP_alongshore/' MOPname '_nowcast.nc'];
    disp('Using nowcast')
    end
elseif date1<tSwitch && date2>tSwitch
    url = ['http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/model/MOP_alongshore/' MOPname '_hindcast.nc'];
    disp('NOTE! Using hindcast dates only')

end

time=ncread(url,'waveTime'); % downloading time from server
time = double(time);
T = datetime(time,'ConvertFrom','posixtime');
disp(['Data is available for ' MOPname  ' from ' datestr(T(1)) ' to ' datestr(T(end))])

itime = find(T>=date1 & T<=date2);
MOP.time = T(itime);
MOP.spec1D=permute(ncread(url,'waveEnergyDensity',[1,min(itime)],[inf,length(itime)]),[2, 1]); % arraging with hours in first value
MOP.frequency = double(ncread(url,'waveFrequency'));
MOP.fbw =  double(ncread(url,'waveBandwidth'));
MOP.Hs = ncread(url,'waveHs',min(itime),length(itime));
MOP.fp = 1./ncread(url,'waveTp',min(itime),length(itime));
MOP.fm = 1./ncread(url,'waveTa',min(itime),length(itime));
MOP.Dp = ncread(url,'waveDp',min(itime),length(itime));

%% get meta data (dimensions of 1, fixed in time)
disp('obtaining meta data')
MOP.depth = ncread(url,'metaWaterDepth');
MOP.lat = ncread(url,'metaLatitude');
MOP.lon = ncread(url,'metaLongitude');
MOP.shorenormal = ncread(url,'metaShoreNormal');
MOP.name = ncread(url,'metaSiteLabel');
MOP.fbounds = ncread(url,'waveFrequencyBounds');

%% QC flags
MOP.flag1 = ncread(url,'waveFlagPrimary',min(itime),length(itime));
MOP.flag2 = ncread(url,'waveFlagSecondary',min(itime),length(itime));

%% directional spreading (2 dimensions)
% Circular moments are calculated from the four lowest fourier coefficients
% a1 b1 a2 b2, following Kuik et al. 1988
disp('Getting directional info with a1 b1 etc.')
a1 = ncread(url,'waveA1Value',[1,min(itime)],[inf,length(itime)]);
b1 = ncread(url,'waveB1Value',[1,min(itime)],[inf,length(itime)]);
a2 = ncread(url,'waveA2Value',[1,min(itime)],[inf,length(itime)]);
b2 = ncread(url,'waveB2Value',[1,min(itime)],[inf,length(itime)]);

m1 = sqrt(a1.^2 + b1.^2);
spread1 = (180/pi) .* ( sqrt( 2 .* ( 1 - m1) ) ); % Eq. 24 <--closest to line equivalent of std
MOP.dirspread1 = nanmean(spread1(1:18,:));

dir2 = atan2(b2,a2)/2; % in radians, Eq. 12
m2 = a2 .* cos(2*dir2) + b2 .* sin(2*dir2); % Eq.20
spread2 = (180/pi) .* sqrt(abs( 0.5 .* ( 1 - m2 ) ) ); % Eq.27 <--less sensitive to reflection
MOP.dirspread2 = nanmean(spread2(1:18,:));
MOP.a1 = a1';
MOP.b1 = b1';
MOP.a2 = a2';
MOP.b2 = b2';

%% radiation stress
MOP.Sxy = ncread(url,'waveSxy',min(itime),length(itime));
MOP.Sxx = ncread(url,'waveSxx',min(itime),length(itime));

end