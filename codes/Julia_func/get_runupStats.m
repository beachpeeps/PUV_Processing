function [Spec,Info,Bulk,Tseries,Sur] = get_runupStats(modDir,inputDate,varargin)
% [Spec,Info,Bulk,Tseries,Sur] = get_runupStats(modDir,inputDate,varargin)
% This function processes the SWASH model data (grd.mat and sur.mat) to get
% the runup line, as defined by Stockdon et al, 2006. 
%
% Input 
% moddir:       file location where sur.mat and grd.mat are located
% inputDate:    not sure why this is needed except for axis labels
% varargin:     additional options for computing the output
%
% Output 
% Spec:     frequency (f), spectra (S), spectra at the offshore
%           boundary (Sboundary), confidence intervals (Slo,Sup), 
%           degrees of freedom (dof);
% Info:     frequency sampling (Hz), runup threshold (threshold),
%           thinLayer***, timestamp in datetime format (datahour), data
%           directory (grdDir, surDir)
% Bulk:     swashparams and confidence intervals, arranged as
%           {Sig,Sinc,eta}, foreshore beach slope (beta), 
%           offshore significant wave height (Ho) 
% Tseries:  timeseries of the runupline (T, Zrunup, Xrunup, idxrunup)
% Sur:      xshore distance (x), mean water level (etabar)
%
% TODO: Is the algorithm finding the correct runupline, especially in
% coarse resolution? Could probably use a little help here.

% default options
options.threshold = 0.05;
options.windowlength = 5; % in minutes, for spectra windows
options.dx = 0.25; % for gridding the runup onto something finer
options.g = 9.81;
options.Tm = 10;
options.wlev = 0;
options.totaltime = 50; % amount to process in minutes, 
options.BoundInd = 1;

options = parseOptions( options , varargin );
BoundInd = options.BoundInd;
%%
% Load the free surface
sur = swash_loadMatFile( [modDir 'sur.mat']);

%Load the grid file (contains Depth and X)
grd  = swash_loadMatFile([modDir 'grd.mat']);

eta = squeeze(sur.Watlev);
Sur.x = grd.Xp;
Sur.etabar = mean(eta,2)+options.wlev;

Wlev = bsxfun(@plus,grd.Botlev,eta');

[nt,~] = size(Wlev);
% put into finer mesh
xx = grd.Xp(1):options.dx:grd.Xp(end); xx = xx(:);
%%
RunupImage = ones(1,nt);
idxrunup = ones(1,nt);
Zrunup = nan(1,nt);
flag = zeros(1,nt);


T = sur.time(1:end); % in seconds

botlevinterp = interp1(grd.Xp,grd.Botlev(1,:),xx);

L2 = [-1000 100; options.threshold options.threshold];

% thinLayer for high slopes when SWASH doesn't do well
    percentData = sum(~isnan(Wlev))./length(sur.time);
    ind80percent = find(percentData>0.8);
    meanRunupShape = nanmean(Wlev);
    thinLayer = meanRunupShape(ind80percent(end));

for i=1:nt
    wlevtemp = Wlev(i,:)';%-thinLayer; % without thinLayer seems to get runup properly now
    
    L1 = [grd.Xp;wlevtemp'];
    runupline = InterX(L1,L2);
    
    if ~isempty(runupline)
        RunupImage(i) = runupline(1);
        idxrunup(i) = find(xx==round(runupline(1)*4)/4);
    else
        % add 0 to wlevtemp for runup intersection ** NOT SURE IF THIS IS
        % BEST PRACTICE **
        % add thin layer
        wlevtemp(find(isnan(wlevtemp),1)) = 0;
        L1 = [grd.Xp;wlevtemp'];
        flag(i) = 1;
        runupline = InterX(L1,L2);
        RunupImage(i) = runupline(1);
        idxrunup(i) = find(xx==round(runupline(1)*4)/4);
    end
    clear wlevtemp L1 runupline
end
%%
Xrunup = RunupImage;
for i=1:nt
    %     XrunupS(i) = xx(RunupImage(i));
    Zrunup(i) = -botlevinterp(idxrunup(i))+options.threshold;
end

Xrunup(RunupImage==1) = nan;
Zrunup(RunupImage==1) = nan;

%%

dt = sur.time(2)-sur.time(1);
dt = floor(dt*10)./10; % want accuracy only to 0.1
%take 1 50 minute chunk for proper spectra comparison
totalSec = floor(options.totaltime*60/dt);

% truncate all records to 50 minutes
Zrunup = Zrunup(end-totalSec-1:end);
Xrunup = Xrunup(end-totalSec-1:end);
idxrunup = idxrunup(end-totalSec-1:end);
T = T(end-totalSec-1:end);

% Get spectra
nfft = options.windowlength*60/dt; % 5 minute chunks
[f, S, Slo, Sup,~,dof] = get_spectrum(detrend(Zrunup(:)), nfft, 1/dt, 0.05);

% find freqs in inc and ig
nIG = find(f>0.004 & f<0.04);
nINC = find(f>=0.04 & f<0.25); % <--careful, do we want to include 0.04 with SS or IG??
df = f(2)-f(1);

%
eta = nanmean(Zrunup(:));
%% get beta
stdEta = std(Zrunup);

maxEta = eta+2*stdEta;
minEta = eta-2*stdEta;

botRange = find(-grd.Botlev>=minEta & -grd.Botlev<=maxEta);

fitvars = polyfit(grd.Xp(botRange), -grd.Botlev(botRange), 1);
beta = fitvars(1);

%% get Ho
[f, Sboundary, ~, ~,~,~] = get_spectrum(detrend(Wlev(:,BoundInd)), nfft, 1/dt, 0.05);
h = grd.Botlev(BoundInd)+Sur.etabar(BoundInd);
k = get_wavenumber(2*pi*f(:),h);
Cg = get_cg(k(:),h);
Cg20 = get_cg(k(:),20);
revshoal = ones(size(f));
revshoal20 = ones(size(f));
ind = find(f<=0.25); % only want to reverse shoal the swell frequencies
revshoal(ind) = 2*Cg(ind)./sqrt(options.g./k(ind));
revshoal20(ind) = Cg(ind)./Cg20(ind);

Soffshore = Sboundary.*revshoal(:);
Soffshore20 = Sboundary.*revshoal20(:);

Hboundary = 4*sqrt(nansum(Sboundary(nINC) * df ) );
Ho_spec = 4*sqrt(nansum(Soffshore(nINC)* df ) );
H20 = 4*sqrt(nansum(Soffshore20(nINC)* df ) );

% Stockdon style, not from spectrum
Ho = reverseshoal(Hboundary,h,1/options.Tm);
g = 9.8;
Lo = options.Tm.^2*g/(2*pi);

% get HoLo_centroid
fswell = [0.04 0.25];
iswell = find(f>=fswell(1) & f<=fswell(2)); % determine swell freqs
fcentroid_swell = sum ( f(iswell).*Soffshore(iswell)' ) ./ sum ( Soffshore(iswell) )' ;
LoC = 9.8./(2*pi.*fcentroid_swell.^2)';



%%
% Sig = 4*sqrt(sum(SLidar(nIG)*df));
% Sinc = 4*sqrt(sum(SLidar(nINC)*df));
% S = sqrt(Sinc.^2 + Sig.^2);

[Sig, Sig_lo, Sig_up, ~ ]= getSWHebounds( S(nIG), dof, 0.95, df );
[Sinc, Sinc_lo, Sinc_up, ~ ]= getSWHebounds( S(nINC), dof, 0.95, df );
Eig = sum(S(nIG).*df);

swashparams = [Sig Sinc eta];
swashparamsLO = [Sig_lo Sinc_lo eta];
swashparamsUP = [Sig_up Sinc_up eta];

Spec.f = f;
Spec.S = S;
Spec.Sboundary = Sboundary;
Spec.Slo = Slo;
Spec.Sup = Sup;
Spec.dof = dof;
Info.Hz = 1/dt;
Info.threshold = options.threshold;
Info.thinLayer = thinLayer;
Info.datahour = datetime(inputDate,'ConvertFrom','datenum');
Info.grdDir =[modDir 'grd.mat'];
Info.surDir = [modDir 'sur.mat'];
Info.waterlevel = options.wlev;
Bulk.swashparams = swashparams;
Bulk.swashparamsLO = swashparamsLO;
Bulk.swashparamsUP = swashparamsUP;
Bulk.swashParamsNames = {'Sig','Sinc','eta'};
Bulk.Eig = Eig;
Bulk.beta = beta;
Bulk.Ho = Ho;
Bulk.Lo = Lo;
Bulk.LoC = LoC;
Bulk.Ho_spec = Ho_spec;
Tseries.T = T; % in seconds
Tseries.Zrunup = Zrunup;
Tseries.Xrunup = Xrunup;
Tseries.idxrunup = idxrunup;






%%
    function options = parseOptions( defaults , cellString )
        
        %INPUT PARSING
        p = inputParser;
        p.KeepUnmatched = true;
        
        names = fieldnames( defaults );
        for ii = 1 : length(names)
            %
            addOptional( p , names{ii}, defaults.(names{ii}));
            %
        end
        %
        parse( p , cellString{:} );
        options = p.Results;
        %
    end

end
