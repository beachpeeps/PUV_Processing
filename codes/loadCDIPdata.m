%% Load data from CDIP website

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

load('buoylist') % fulll list of all CDIP buoy's as of 07/22/19
[indx,tf] = listdlg('ListString',buoylist, 'SelectionMode','single', 'InitialValue',[1], 'Name', 'What station do you want to use?');
stn = char(buoylist(indx));
stn = string(stn(end-2:end))

st = str2double(stn);
%CONNECT TO THREDDS SERVER AND OPEN NETCDF FILE
if st == 073
    urlbase = 'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/xarchive/';
elseif st == 239 || st == 238 || st ==245 || st ==240
    fprintf('No historic archived data avaliable')
    return
else 
    urlbase = 'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/archive/';
end
p1 = 'p1/';
urlend ='p1_historic.nc';
dsurl = char(strcat(urlbase,stn,p1,stn,urlend));
ds= ncdataset(dsurl);

%GET TIMES OF DATA
timevar = ds.data('waveTime');
timeconvert = ds.time('waveTime',timevar); % Convert UNIX timestamps to Matlab serial units
timeseriesbuoy = datetime(timeconvert, 'ConvertFrom', 'datenum');

Hsbuoy = ds.data('waveHs');
Ed = ds.data('waveEnergyDensity');
Fq = ds.data('waveFrequency');
Bw = ds.data('waveBandwidth')';
Dp = ds.data('waveDp');
Tp = ds.data('waveTp');

%% finding peak buoy Hs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[peaks, loc] = findpeaks(Hsbuoy);
maxpeaksbuoy = maxk(peaks, 6)
for i = 1:size(maxpeaksbuoy)
    maxlocbuoy(i) = loc(find(maxpeaksbuoy(i)== peaks))
end

for i = 1:6
    for j = i+1:6
        if abs(maxlocbuoy(i)-maxlocbuoy(j)) < 100
           maxlocbuoy(j) = 0;
        else 
        end
    end
end
maxpeaksbuoy = maxpeaksbuoy(maxlocbuoy~=0)
maxlocbuoy = maxlocbuoy(maxlocbuoy~=0)
timeseriesbuoy(maxlocbuoy) %gives times of maximum

