%% SAVEING RAW DATA TO NETCDF - LOOPS THROUGH ALL 'SETS' OF DATA FOUND IN FOLDER
%
%
%
% Author: 
% Athina Lange, SIO July 2021

for ii = 1:length(filenames)
    eval(['DAT = DAT' char(string(ii)) ';']);
    eval(['SEN = SEN' char(string(ii)) ';']);
    
    newfilename = [filename '_raw_' char(string(ii)) '.nc'];
    
    coordsys = 'XYZ'


    % set up netcdf

    ncid = netcdf.create(newfilename, 'NC_WRITE');
    dimidlat = netcdf.defDim(ncid,'latitude',1);
    dimidlon = netcdf.defDim(ncid,'longitude',1);
    dimidrot_sensor = netcdf.defDim(ncid,'sensor_rotation',1);
    dimidfs = netcdf.defDim(ncid,'fs',1);
    dimidclock = netcdf.defDim(ncid,'clockdrift',1);
    dimidcoordsys = netcdf.defDim(ncid,'coordsys',length(coordsys));

    dimidDAT = netcdf.defDim(ncid, 'dat_size',size(DAT,2));
    dimidSEN = netcdf.defDim(ncid, 'sen_size',size(SEN,2));

    dimid_DATtime = netcdf.defDim(ncid,'time_dat',size(DAT,1));
    dimid_SENtime = netcdf.defDim(ncid,'time_sen',size(SEN,1));

    
    latitude_ID=netcdf.defVar(ncid,'latitude','double',[dimidlat]);
    longitude_ID=netcdf.defVar(ncid,'longitude','double',[dimidlon]);
    rot_sensor_ID=netcdf.defVar(ncid,'sensor_rotation','double',[dimidrot_sensor]);
    fs_ID=netcdf.defVar(ncid,'fs','double',[dimidfs]);
    clock_ID=netcdf.defVar(ncid,'clockdrift','double',[dimidclock]);
    coordsys_ID=netcdf.defVar(ncid,'coordsys','char',[dimidcoordsys]);

    DAT_ID = netcdf.defVar(ncid,'DAT','double', [dimid_DATtime dimidDAT]);
    SEN_ID = netcdf.defVar(ncid,'SEN','double', [dimid_SENtime dimidSEN]);

    netcdf.endDef(ncid);
    %%
    netcdf.putVar(ncid,latitude_ID,LAT);
    netcdf.putVar(ncid,longitude_ID,LON);
    netcdf.putVar(ncid,rot_sensor_ID,rot_angle);
    netcdf.putVar(ncid,fs_ID,fs);
    netcdf.putVar(ncid,clock_ID,clockdrift);
    netcdf.putVar(ncid,coordsys_ID,coordsys);

    netcdf.putVar(ncid, DAT_ID, DAT);
    netcdf.putVar(ncid, SEN_ID, SEN);
    netcdf.close(ncid);
    %%

    ncwriteatt(newfilename, 'longitude', 'units', 'degrees_east');
    ncwriteatt(newfilename, 'latitude', 'units', 'degrees_north');

    ncwriteatt(newfilename, 'sensor_rotation', 'long_name', 'Heading of Beam1 (Compass)');
    ncwriteatt(newfilename, 'fs', 'long_name', 'Sampling Frequency');
    ncwriteatt(newfilename, 'fs', 'units', 'Hz');
    
    ncwriteatt(newfilename, 'clockdrift', 'long_name', 'Clock drift from excel');
    ncwriteatt(newfilename, 'clockdrift', 'units', 's');



    ncwriteatt(newfilename, 'DAT', 'long_name', 'contents of .dat file');
    ncwriteatt(newfilename, 'DAT', 'Column 1', '1 Burst counter');
    ncwriteatt(newfilename, 'DAT', 'Column 2', '2 Ensemble counter');
    ncwriteatt(newfilename, 'DAT', 'Column 3', '3 Velocity (X) m/s');
    ncwriteatt(newfilename, 'DAT', 'Column 4', '4 Velocity (Y) m/s');
    ncwriteatt(newfilename, 'DAT', 'Column 5', '5 Velocity (Z) m/s');
    ncwriteatt(newfilename, 'DAT', 'Column 6', '6 Amplitude (Beam 1) counts');
    ncwriteatt(newfilename, 'DAT', 'Column 7', '7 Amplitude (Beam 2) counts');
    ncwriteatt(newfilename, 'DAT', 'Column 8', '8 Amplitude (Beam 3) counts');
    ncwriteatt(newfilename, 'DAT', 'Column 9', '9 SNR (Beam 1) dB');
    ncwriteatt(newfilename, 'DAT', 'Column 10', '10 SNR (Beam 2) dB');
    ncwriteatt(newfilename, 'DAT', 'Column 11', '11 SNR (Beam 3) dB');
    ncwriteatt(newfilename, 'DAT', 'Column 12', '12 Correlation (Beam 1) %');
    ncwriteatt(newfilename, 'DAT', 'Column 13', '13 Correlation (Beam 2) %');
    ncwriteatt(newfilename, 'DAT', 'Column 14', '14 Correlation (Beam 3) %');
    ncwriteatt(newfilename, 'DAT', 'Column 15', '15 Pressure dBar');
    ncwriteatt(newfilename, 'DAT', 'Column 16', '16 Analog input 1');
    ncwriteatt(newfilename, 'DAT', 'Column 17', '17 Analog input 2');
    ncwriteatt(newfilename, 'DAT', 'Column 18', '18 Checksum (1=failed)');
    
    ncwriteatt(newfilename, 'SEN', 'long_name', 'contents of .sen file');
    ncwriteatt(newfilename, 'SEN', 'Column 1', '1 Month (1-12)');
    ncwriteatt(newfilename, 'SEN', 'Column 2', '2 Day (1-31)');
    ncwriteatt(newfilename, 'SEN', 'Column 3', '3 Year');
    ncwriteatt(newfilename, 'SEN', 'Column 4', '4 Hour (0-23)');
    ncwriteatt(newfilename, 'SEN', 'Column 5', '5 Minute (0-59)');
    ncwriteatt(newfilename, 'SEN', 'Column 6', '6 Second (0-59)');
    ncwriteatt(newfilename, 'SEN', 'Column 7', '7 Error code');
    ncwriteatt(newfilename, 'SEN', 'Column 8', '8 Status Code');
    ncwriteatt(newfilename, 'SEN', 'Column 9', '9 Battery voltage');
    ncwriteatt(newfilename, 'SEN', 'Column 10', '10 Soundspeed m/s');
    ncwriteatt(newfilename, 'SEN', 'Column 11', '11 Heading deg');
    ncwriteatt(newfilename, 'SEN', 'Column 12', '12 Pitch deg');
    ncwriteatt(newfilename, 'SEN', 'Column 13', '13 Roll deg');
    ncwriteatt(newfilename, 'SEN', 'Column 14', '14 Temperature degC');
    ncwriteatt(newfilename, 'SEN', 'Column 15', '15 Analog input');
    ncwriteatt(newfilename, 'SEN', 'Column 16', '16 Checksum (1=failed)');
   
    ncdisp(newfilename)
end
