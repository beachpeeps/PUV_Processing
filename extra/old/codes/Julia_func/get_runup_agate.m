function [Spec,Info,Bulk,Tseries] = get_runup_agate(folderdir,filedir,filename,casenum,threshold)
%%
if exist('threshold','var') == 0
    threshold = 0.1;
end

[res,datahour,~,wlev,Wlev,X,M,N, bulkinfo]...
    = get_matfile(folderdir,filedir,casenum);

%%
U = squeeze(res.Vksi(1,:,:))';
Wlev = bsxfun(@plus,res.Botlev,squeeze(res.Watlev(1,:,:))');
% check = squeeze(res.Watlev(1,:,:))';
% Wlev1 = squeeze(res.Watlev(1,:,:))';
%% THIS SECTION WAS A MISTAKE, 
% KEEP FOR POSTERITY AND REMINDER never to do this again
% if reflection
%     fmin = 0.001; fmax=0.04;
%     H = mean(Wlev1);
%     
%     for i=1:M
%         [Wlev(:,i),~] = RemoveReflection( check(:,i), U(:,i), H(i),  0.5,'method','sher','flim', [fmin,fmax],'g',9.81,'trend',true,'depav',true);
%     end
%     Wlev2 = bsxfun(@plus,res.Botlev,Wlev);
%     
% elseif ~reflection
%     Wlev = Wlev1;
% end
%%
% Wlev = Wlev';
% [M,N] = size(Wlev);
% X = res.Xp(1,:); 
X= X(:);
xx = [-1400:0.25:-100]; xx = xx(:);
%%

%ensure N is even
if floor(N/2) ~= N/2
    N = N-1;
end


% Wlev = Wlev+res.Botlev(1,:)'*ones(1,N);
% Wlev = bsxfun(@plus,res.Botlev,Wlev);
Wlev = Wlev';
Wlev = Wlev(:,1:N);

%%
RunupImage = ones(1,N);
idxrunup = ones(1,N);
Xrunup = nan(1,N);
Zrunup = nan(1,N);
T = datahour+ res.time(1:end)./3600/24;
T = T(1:N);
botlevinterp = interp1(X,res.Botlev(1,:),xx);
% wlevinterp = nan(length(xx),N);
%%
L2 = [-1000 -100; threshold threshold];

for i=1:N
    wlevtemp = Wlev(:,i);
%     wlevinterp = interp1(X,wlevtemp,xx);

    L1 = [X';wlevtemp'];
    runupline = InterX(L1,L2);
    runupline = fliplr(runupline);

%     wlevinterp(:,i) = interp1(X,wlevtemp,xx);
%     runupline = find(wlevinterp(:,i)<=0.1,1);
    
    if ~isempty(runupline)
%         RunupImage(i) = runupline(1);
        RunupImage(i) = round(runupline(1)*4)/4;
        idxrunup(i) = find(xx==RunupImage(i));

    end
end
%%
XrunupS = RunupImage;
for i=1:N
%     XrunupS(i) = xx(RunupImage(i));
    ZrunupS(i) = wlev-botlevinterp(idxrunup(i));
end

XrunupS(RunupImage==1) = nan;
ZrunupS(RunupImage==1) = nan;

%%

Hz = res.time(2)-res.time(1);
Hz = floor(Hz*10)./10;
%take 1 50 minute chunk for proper spectra comparison
totalSec = floor(50*60/Hz);


ZrunupS = ZrunupS(end-totalSec-1:end);
XrunupS = XrunupS(end-totalSec-1:end);
idxrunup = idxrunup(end-totalSec-1:end);

% ZrunupS = ZrunupS(1:totalSec);

T = T(1:length(ZrunupS));
znonan = inpaint_nans(ZrunupS,4);
% zsmooth = smooth(znonan,5);
zsmooth = znonan';


% pp = 12;
% [PSI,~] = sleptap(length(znonan),pp);
% [fLidar,SLidar] = mspec(Hz,detrend(zsmooth),PSI);
nfft = 10*60;
[fLidar, SLidar, SLlo, SLup,~,dof] = get_spectrum(detrend(zsmooth), nfft, 1/Hz, 0.05);

% For computing Stockdon eq. 7: Swash = C*sqrt(S_inc^2 + S_ig^2)

% find freqs in inc and ig
nIG = find(fLidar>0.004 & fLidar<0.04);
nINC = find(fLidar>0.04 & fLidar<0.25);
df = fLidar(2)-fLidar(1);

%
eta = nanmean(ZrunupS(:));
% Sig = 4*sqrt(sum(SLidar(nIG)*df));
% Sinc = 4*sqrt(sum(SLidar(nINC)*df));
% S = sqrt(Sinc.^2 + Sig.^2);

[Sig, Sig_lo, Sig_up, ~ ]= getSWHebounds( SLidar(nIG), dof, 0.95, df );
[Sinc, Sinc_lo, Sinc_up, ~ ]= getSWHebounds( SLidar(nINC), dof, 0.95, df );


swashparams = [Sig Sinc eta];
swashparamsLO = [Sig_lo Sinc_lo eta];
swashparamsUP = [Sig_up Sinc_up eta];

ZrunupS = zsmooth;

Spec.f = fLidar;
Spec.S = SLidar;
Spec.Slo = SLlo;
Spec.Sup = SLup;
Spec.dof = dof;
Info.datahour = datahour;
Info.Hz = Hz;
Info.bulkinfo = bulkinfo;
Info.wlev = wlev;
Bulk.swashparams = swashparams;
Bulk.swashparamsLO = swashparamsLO;
Bulk.swashparamsUP = swashparamsUP;
Tseries.T = T;
Tseries.Zrunup = ZrunupS;
Tseries.Xrunup = XrunupS;
Tseries.idxrunup = idxrunup;


