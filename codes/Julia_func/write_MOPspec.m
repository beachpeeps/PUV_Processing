function [Hs, dir, fspread, sqrtHoLo, boundcondfile] = write_MOPspec(MOPname,t,outdir,BC)
% write_MOPspec(MOPname,t)
% INPUT:
%   MOPname: example: 'D0045'
%   t: datetime vector
%   outdir: directory where file will be written
%   Eig: 'Eig=0' for zero IG energy at boundary
%
% OUTPUT:
%   1D spectrum file for SWASH input, named by t input, ISO 8601 compliant
%   convention, yyyymmddTHHMMSS. Time is same as MOP time, in UTC.
%   Hs: significant wave height 
%   dir: peak wave direction
%   fspread: frequency spread
%   sqrtHoLo: sqrt(deepwater H * wavelength)


%% read input
% MOP = read_MOPline(MOPname,t,t);
load('~/Documents/NeoStockdon/data/binnedMOPS.mat','MOP')
tind = find(t==MOP.time);

%% prep data
% artifically decrease df so no waves are repeated
% TODO: is this acceptable??
cycle = 60*60; % set non-repeating cycle to 60 minutes, in seconds
df = 1/cycle;
freq = MOP.frequency;
dfreq = double(MOP.fbounds(2,:)-MOP.fbounds(1,:))'; % f bounds are not constant!!

f = freq(1) : df : freq(end); 

S = interp1(freq,MOP.spec1D(tind,:),f); % interpolate coarse E(f) spec to finer scale

%% add 0 to IG band for Eig = 0
switch BC
    case 'Eig=0'
        ff = 0.004:df:f(1)-df;
        SS = zeros(size(ff));

        f = [ff f];
        S = [SS S];
end
%% calculate output
Hs = MOP.Hs(tind);
dir = MOP.Dp(tind);
fswell = [0.04 0.25];
iswell = find(freq>=fswell(1) & freq<=fswell(2)); % determine swell freqs

% centroid frequency
fc = sum ( freq(iswell).* MOP.spec1D(tind,iswell)' ) ./ sum ( MOP.spec1D(tind,iswell),2 )' ;

% frequency spread: TODO: coding on this is insanely hard to read 
fspread = sqrt(sum ( (freq(iswell)-fc).^2.* MOP.spec1D(tind,iswell)'.*dfreq(iswell) ) ./ sum ( MOP.spec1D(tind,iswell).*dfreq(iswell)',2 )' );

% sqrt( Ho Lo)
Tp = 1./MOP.fp(tind);
Lo = 9.8/(2*pi)*Tp.^2;
sqrtHoLo = sqrt(MOP.Hs(tind).*Lo);

%% write to file
specname = datestr(MOP.time(tind),'yyyymmddTHH'); % naming convention 'yyyymmddTHHMMSS'
boundcondfile = [specname '.spec'];
fid   = fopen([outdir boundcondfile],'w'); % open file to write
fprintf(fid,'%s\n','SPEC1D'); % first line is header 'SPEC1D'


for i=1:length(f)
    fprintf(fid,'%10.4E %10.4E\n',f(i), S(i)); % write spectrum
end
fclose('all');

disp(['Created ' outdir boundcondfile])
end



%
%
