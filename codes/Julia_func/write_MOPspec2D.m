function [Hs, dir, fspread, sqrtHoLo, boundcondfile] = write_MOPspec2D(xpinp,MOPname,t,outfiledir)
%% read input
MOP = read_MOPline(MOPname,t,t);
Loc = [xpinp, 0]; % [x,y], puts location at the vertical middle and west side boundary

Dir = 1:360;
Fname = 'wave_forcing_2D.bnd'; %TODO: fix for variable name
%%
freq = MOP.frequency;

lfreq = length(freq);
ldir = length(Dir);

%% The following lines use the MEM estimator to determine a1 b1 a2 b2. 
mem = zeros(lfreq,ldir);
for fi=1:lfreq  
    mem(fi,:) = get_mem(MOP.a1(fi),MOP.b1(fi),MOP.a2(fi),MOP.b2(fi));
end

%% need to rotate angle into shorenormal coordinates
% is this right??? TODO: fix this!!
% rotateBy = round(270-MOP.shorenormal);
% cutAt = 270+rotateBy;

rotateBy = round(270-MOP.shorenormal);
cutAt = 90+rotateBy;
mem = mem(:,[cutAt:end 1:cutAt-1]);

%% multiply to get E(f,theta)
Sf      = repmat(MOP.spec1D',[1 360]);
Sfth    = mem.*Sf;

% interpolate to finer df mesh to avoid cycling issues
cycle = 60*60; % set non-repeating cycle to 60 minutes, in seconds
df = 1/cycle;
f = freq(1) : df : freq(end); 
lf = length(f);

Sold = griddedInterpolant({freq,Dir},Sfth);
Snew = Sold({f,Dir});

%%
Spec = zeros(1,lf,ldir); %Spec(NLOC,NFREQ,NDIR)
Spec(1,:,:) = Snew;

%% [YYYY,MO,DD,HH,MM,SS] = datevec(t(i));
boundcondfile = 'wave_forcing_2D.bnd';
createswash2Dspec(Loc,f,Dir,Spec,outfiledir)

%% calculate output
Hs = MOP.Hs;
dir = MOP.Dp;
fswell = [0.04 0.25];
iswell = find(freq>=fswell(1) & freq<=fswell(2)); % determine swell freqs
dfreq = double(MOP.fbounds(2,:)-MOP.fbounds(1,:))'; % f bounds are not constant!!

% centroid frequency
fc = sum ( freq(iswell).* MOP.spec1D(iswell)' ) ./ sum ( MOP.spec1D(iswell),2 )' ;

% frequency spread: TODO: coding on this is insanely hard to read 
fspread = sqrt(sum ( (freq(iswell)-fc).^2.* MOP.spec1D(iswell)'.*dfreq(iswell) ) ./ sum ( MOP.spec1D(iswell).*dfreq(iswell)',2 )' );

% sqrt( Ho Lo)
Tp = 1./MOP.fp;
Lo = 9.8/(2*pi)*Tp.^2;
sqrtHoLo = sqrt(MOP.Hs.*Lo);