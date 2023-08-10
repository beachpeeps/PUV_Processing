%% Load data from MOP website

% allows you to choose which buoy's data you want to look at
% (from list of all avaliable CDIP buoys (list updated 07/22/19))
% pull historical data as a netCDF file (ds)
% make sure the nctoolbox is installed in the MATLAB path
% CAUTION: different buoys have different url's to their historical data ->
%       code works for stn 121 (Ipan), but not all station's have been checked
% creates a datetime vector of all the instances that there is avaliable
% data from the buoy (timeseries)

%make sure the nctoolbox is installed in the MATLAB path
%may need to change url depending on station

%clear all
setup_nctoolbox

load('moplist') % fulll list of all CDIP buoy's as of 07/22/19
[indx,tf] = listdlg('ListString',moplist, 'SelectionMode','single', 'InitialValue',[1], 'Name', 'What MOP do you want to use?');
stn = char(moplist(indx));
stn = string(stn(end-4:end))

st = str2double(stn);
%CONNECT TO THREDDS SERVER AND OPEN NETCDF FILE
%if st == 073
%    urlbase = 'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/xarchive/';
%elseif st == 239 || st == 238 || st ==245 || st ==240
%    fprintf('No historic archived data avaliable')
%    return
%else 
    urlbase = 'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/model/MOP_alongshore/';
%end
urlend ='_nowcast.nc';
ncfile = strcat(urlbase,stn,urlend);


%GET TIMES OF DATA
timeseries_mop = datetime(double(ncread(ncfile,'waveTime'))/86400 + datenum(1970,1,1), 'ConvertFrom', 'datenum');

Hs_mop = ncread(ncfile, 'waveHs');
Ed_mop = ncread(ncfile, 'waveEnergyDensity');
Fq_mop = ncread(ncfile, 'waveFrequency');
Bw_mop = ncread(ncfile, 'waveBandwidth')';
Dp_mop = ncread(ncfile, 'waveDp');
Tp_mop = ncread(ncfile, 'waveTp');
Sxx_mop = ncread(ncfile, 'waveSxx');
Sxy_mop = ncread(ncfile, 'waveSxy');
Dirmean_mop = ncread(ncfile, 'waveMeanDirection');
a1_mop = ncread(ncfile, 'waveA1Value');
b1_mop = ncread(ncfile, 'waveB1Value');
a2_mop = ncread(ncfile, 'waveA2Value');
b2_mop = ncread(ncfile, 'waveB2Value');
snorm=ncread(ncfile,'metaShoreNormal'); % mop shore normal

