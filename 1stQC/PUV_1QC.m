%% Import raw PUV data and do first QC
%
%
%
% CAUTION: work still TBD for input beginning & inspection times -
%           currently values hard coded
%
%
%% Author: 
% Athina Lange, SIO July 2021

tic
% input directory where data is - one deloyment per folder
%directory = '/Volumes/Lange_SIO/Ch2_IG/DATA/By instrument/NortekVector/2019-2020_TorreyPinesMOP582_10meter'
directory = '/Users/athinalange/Desktop/IG_BC_remote_work/2019-2020_TorreyPinesMOP582_10meter/DATA'
directory = '~/Desktop/check_Athina_code/2019-2020_TorreyPinesMOP582_10meter'

filename = 'Torrey_1920_582_10'

%% From fieldwork excel sheet:
% input lat/lon of sensor
LATLON = [32.926405	-117.265776]; % Pull info from deployment files

% input compass heading of beam 1 of ADV
rot_angle = [57.3]; % rotation to magnetic North (from excel sheets)

% input clock drift between beginning and end of deployment
clockdrift = 7; % sec


%% run raw_processing codes
[PUV] = PUV_raw_process(directory, filename, LATLON, rot_angle, clockdrift)

toc