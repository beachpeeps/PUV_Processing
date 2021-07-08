%% Extract data from .dat files
tic
% input directory where data is - one deloyment per folder
%directory = '/Volumes/Lange_SIO/Ch2_IG/DATA/By instrument/NortekVector/2019-2020_TorreyPinesMOP582_10meter'
directory = '/Users/athinalange/Desktop/IG_BC_remote_work/2019-2020_TorreyPinesMOP582_10meter'
filename = 'Torrey_1920_582_10'

% input lat/lon of sensor
LATLON = [32.926405	-117.265776]; % Pull info from deployment files
% input compass heading of beam 1 of ADV
rot_angle = [57.3]; % rotation to magnetic North (from excel sheets)
% input clock drift between beginning and end of deployment
clockdrift = 7; % sec

[PUV] = PUV_raw_process(directory, filename, LATLON, rot_angle, clockdrift)

toc
