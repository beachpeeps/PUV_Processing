%% SAVE PUV DATA TO NETCDF AFTER INITIAL QC

newfilename = [filename '_processed.nc'];

t = datenum(PUV.time)-datenum(2000,1,1);


% set up netcdf

ncid = netcdf.create(newfilename, 'NC_WRITE');
dimidt = netcdf.defDim(ncid,'time',length(t));
dimidlat = netcdf.defDim(ncid,'latitude',1);
dimidlon = netcdf.defDim(ncid,'longitude',1);
dimidrot_sensor = netcdf.defDim(ncid,'sensor_rotation',1);
dimidrot_mag = netcdf.defDim(ncid,'mag_rotation',1);
dimidfs = netcdf.defDim(ncid,'fs',1);


time_ID=netcdf.defVar(ncid,'time','double',[dimidt]);
latitude_ID=netcdf.defVar(ncid,'latitude','double',[dimidlat]);
longitude_ID=netcdf.defVar(ncid,'longitude','double',[dimidlon]);
rot_sensor_ID=netcdf.defVar(ncid,'sensor_rotation','double',[dimidrot_sensor]);
rot_mag_ID=netcdf.defVar(ncid,'mag_rotation','double',[dimidrot_mag]);
fs_ID=netcdf.defVar(ncid,'fs','double',[dimidfs]);


P_ID = netcdf.defVar(ncid,'P','double',[dimidt]);
U_ID = netcdf.defVar(ncid,'U','double',[dimidt]);
V_ID = netcdf.defVar(ncid,'V','double',[dimidt]);
W_ID = netcdf.defVar(ncid,'W','double',[dimidt]);
U_xy_ID = netcdf.defVar(ncid,'U_xy','double',[dimidt]);
V_xy_ID = netcdf.defVar(ncid,'V_xy','double',[dimidt]);
T_ID = netcdf.defVar(ncid,'T','double',[dimidt]);
netcdf.endDef(ncid);


netcdf.putVar(ncid,time_ID,t);
netcdf.putVar(ncid,latitude_ID,PUV.LATLON(1));
netcdf.putVar(ncid,longitude_ID,PUV.LATLON(2));
netcdf.putVar(ncid,rot_sensor_ID,PUV.rotation.sensor);
netcdf.putVar(ncid,rot_mag_ID,PUV.rotation.mag);
netcdf.putVar(ncid,fs_ID,PUV.fs);

netcdf.putVar(ncid,P_ID,PUV.P);
netcdf.putVar(ncid,U_ID,PUV.BuoyCoord.U);
netcdf.putVar(ncid,V_ID,PUV.BuoyCoord.V);
netcdf.putVar(ncid,W_ID,PUV.BuoyCoord.W);
netcdf.putVar(ncid,U_xy_ID,PUV.InstrCoord.U);
netcdf.putVar(ncid,V_xy_ID,PUV.InstrCoord.V);
netcdf.putVar(ncid,T_ID,PUV.T);
netcdf.close(ncid);

ncwriteatt(newfilename, 'time', 'units', 'days since 2000-01-01 00:00:00');
ncwriteatt(newfilename, 'time', 'long_name', 'julian day (UT)');
ncwriteatt(newfilename, 'time', 'conventions', 'relative julian days with decimal part (as parts of the day )');
ncwriteatt(newfilename, 'time', 'calendar', 'standard');
ncwriteatt(newfilename, 'time', 'notes', 'adjusted for clockdrift');

ncwriteatt(newfilename, 'longitude', 'units', 'degrees_east');
ncwriteatt(newfilename, 'latitude', 'units', 'degrees_north');

ncwriteatt(newfilename, 'sensor_rotation', 'long_name', 'Heading of Beam1 (Compass)');
ncwriteatt(newfilename, 'mag_rotation', 'long_name', 'Magnetic Declination');
ncwriteatt(newfilename, 'mag_rotation', 'model', 'MATLAB igrfmagm - IGRF-13');
ncwriteatt(newfilename, 'fs', 'long_name', 'Sampling Frequency');
ncwriteatt(newfilename, 'fs', 'units', 'Hz');


ncwriteatt(newfilename, 'P', 'units', 'dBar');
ncwriteatt(newfilename, 'P', 'QC', 'Beginning of deployment, Inspection, and Battery depletion times removed; min(Corr)> 70%;');

ncwriteatt(newfilename, 'U', 'units', 'm s-1');
ncwriteatt(newfilename, 'U', 'axis', 'X_rotated');
ncwriteatt(newfilename, 'U', 'positive', 'WEST');
ncwriteatt(newfilename, 'U', 'QC', 'Beginning of deployment, Inspection, and Battery depletion times removed; min(Corr)> 70%;');
ncwriteatt(newfilename, 'V', 'units', 'm s-1');
ncwriteatt(newfilename, 'V', 'axis', 'Y_rotated');
ncwriteatt(newfilename, 'V', 'positive', 'NORTH');
ncwriteatt(newfilename, 'V', 'QC', 'Beginning of deployment, Inspection, and Battery depletion times removed; min(Corr)> 70%;');
ncwriteatt(newfilename, 'W', 'units', 'm s-1');
ncwriteatt(newfilename, 'W', 'axis', 'Z_rotated');
ncwriteatt(newfilename, 'W', 'positive', 'UP');
ncwriteatt(newfilename, 'W', 'QC', 'Beginning of deployment, Inspection, and Battery depletion times removed; min(Corr)> 70%;');

ncwriteatt(newfilename, 'U_xy', 'units', 'm s-1');
ncwriteatt(newfilename, 'U_xy', 'axis', 'X');
ncwriteatt(newfilename, 'U_xy', 'positive', 'Cross-shore normal');
ncwriteatt(newfilename, 'U_xy', 'QC', 'Beginning of deployment, Inspection, and Battery depletion times removed; min(Corr)> 70%;');
ncwriteatt(newfilename, 'V_xy', 'units', 'm s-1');
ncwriteatt(newfilename, 'V_xy', 'axis', 'Y');
ncwriteatt(newfilename, 'V_xy', 'positive', 'Alongshore south');
ncwriteatt(newfilename, 'V_xy', 'QC', 'Beginning of deployment, Inspection, and Battery depletion times removed; min(Corr)> 70%;');

ncwriteatt(newfilename, 'T', 'units', 'C');
ncwriteatt(newfilename, 'T', 'QC', 'Beginning of deployment, Inspection, and Battery depletion times removed; min(Corr)> 70%;');

ncdisp(newfilename)
