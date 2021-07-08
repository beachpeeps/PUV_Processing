# PUV_Processing
CPG PUV Processing Codes

Run PUV_1QC.m

Change data directory and filename, pull lat & lon, heading angle and clockdrift from Excel sheet. 


Will save raw data into netcdf files.

Will save qc'd data into 1 netcdf file & .mat file (as structure).


[PUV.ShorenormalCoord.U, PUV.ShorenormalCoord.V, shorenormal] = rotate_shorenormal(582, PUV.BuoyCoord.U, PUV.BuoyCoord.V)

