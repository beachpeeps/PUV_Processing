function createswash2Dspec(Loc,Freq,Dir,Spec,outfiledir)
%function createswash2Dspec(Freq,Dir,Spec,Fname)
%USAGE: This function takes a directional spectrum and outputs it into a
%boundary forcing file suitable for SWASH.
% 
%INPUT:
%Loc     = Co-ordinates for location at which boundary information is provided (i.e., Loc = [Lon,Lat])
%Freq    = Vector of frequency spacing of the two-dimensional spectrum
%Dir     = Vector of directional spacing of the two-dimensional spectrum
%Spec    = Frequency-directional spectrum (expected units is m^2/Hz/Deg). Also it is assumed that the
%          spectrum is a 2-dimensional variable with a structure Spec(NFREQ,NDIR).
%          NFREQ is the number of frequencies and NDIR is the number of directions
%Fname   = File name for output

% Nirnimesh Kumar and George Voulgaris
% Coastal Processes and Sediment Dynamics Lab
% Dept. of Earth and Ocean Sciences,
% Univ. of South Carolina, Columbia, SC
% 03/06/2013
% Adapted by Julia Fiedler for use in SWASH model 4/19/2019

if (nargin<5 || isempty(outfiledir)==1)
    Fname = 'wave_forcing.bnd';
end

if nargin<4
    error('Function requires location, frequency, direction and spec')
end

Fname = [outfiledir 'wave_forcing_2D.bnd'];

NLOC  = 1               ;     %Number of locations
NFREQ = length(Freq)    ;     %Number of frequencies
NDIR  = length(Dir)     ;     %Number of directions
rho   = 1025            ;     %SWAN allows changing rho value, but deafult is 1025 kg/m^3
g     = 9.81            ;     %Accn. due to gravity
ver   = '999'			;     %SWAN version
proj  = 'NeoStockdon'   ;	  %Project Name
Exval = -99				;	  %Exception value

fid   = fopen(Fname,'w');
fprintf(fid,'%s\n','SWAN   1                                Swan standard spectral file, version');
fprintf(fid,'%s\n',['$   Data produced by SWAN version ',ver]);
fprintf(fid,'%s\n',['$   Project: ',proj,'        ;  run number:']);
fprintf(fid,'%s\n','LOCATIONS                                  locations in x-y space');
fprintf(fid,'%6i',NLOC);
fprintf(fid,'%s\n','                                  number of locations');
for i=1:1:NLOC
    fprintf(fid,'%12.6f %12.6f\n',Loc(i,:));
end
clear i

fprintf(fid,'%s\n','AFREQ                                   absolute frequencies in Hz');
fprintf(fid,'%6i',NFREQ);
fprintf(fid,'%s\n','                                  number of frequencies');
for i=1:1:NFREQ
    fprintf(fid,'%10.4f\n',Freq(i));
end
clear i
fprintf(fid,'%s\n','NDIR                                    spectral nautical directions in degr');
fprintf(fid,'%6i',NDIR);
fprintf(fid,'%s\n','                                  number of directions');
for i=1:1:NDIR
    if Dir(i)>360
        Dir(i)=Dir(i)-360;
    end
    fprintf(fid,'%10.4f\n',Dir(i));
end
clear i
fprintf(fid,'%s\n', 'QUANT');
fprintf(fid,'%s\n', '     1                                  number of quantities in table');
fprintf(fid,'%s\n', 'VaDens                                  variance densities in m2/Hz/degr');
fprintf(fid,'%s\n', 'm2/Hz/degr                             unit');
fprintf(fid,'%14.4e',Exval);
fprintf(fid,'%s\n', '                          exception value');

for q=1:1:NLOC
    Edens = squeeze(Spec(q,:,:));
    [ID]=max(max(Edens));
    if ID>(10^-10)
        FAC=(1.01*ID*10^(-4));
    else
        FAC=10^(-5);
    end
    Edens = Edens/FAC;
    FACS  = find_exp(FAC);
    fprintf(fid,'%s\n','FACTOR');
    fprintf(fid,'%s\n',FACS);
    for r=1:1:NFREQ
        fprintf(fid,'%6d',round(Edens(r,:)));
        fprintf(fid,'\n');
    end
end

fclose(fid);
end

function [N]=find_exp(x);
%
% Convert number x into a string with an exponent format
% xx.xxxxxxxxE-nn

i=-1;
if floor(x)==0
    g=0;
    while g<1
        i=i+1;
        g=floor(rem(x*(10^i),10));
    end
    N = -i;
    G=x*(10^-N);
else
    g=500;
    while (g~=0)
        i=i+1;
        g=floor(x/(10^i));
    end
    N = i-1;
    G=x*(10^-N);
end
N=sprintf('%13.8f%s%+03i',G,'E',N);
end



