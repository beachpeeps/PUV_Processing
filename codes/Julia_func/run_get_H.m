clear
volumedir = '/Volumes/jfiedler/SWASH_runup/';
inFileDir = [volumedir 'mat/bathytest0_005f/'];
outFileDir = [volumedir 'mat/processed/'];
datahour = datenum(2013,10,21,13,0,0);

get_H_bathytests(inFileDir,outFileDir,datahour)