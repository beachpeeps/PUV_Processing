function [freq,Snn] = read_swash_spec(inputFolder)
% % [xgrid, ygrid] = read_swash_coord(inputFile)
% % function to read swash_coord.grd file located in inputFolder
% % This should be universal to any grid but has not been thoroughly tested
% 
% % Julia Fiedler March 2019
% 
% % get input
inputFile = dir([inputFolder '*.spec']);
fID = fopen([inputFolder inputFile.name]);

C = textscan(fID,'%f %f','HeaderLines',2,'Delimiter',' ');

freq = C{1,1};
Snn = C{1,2};

