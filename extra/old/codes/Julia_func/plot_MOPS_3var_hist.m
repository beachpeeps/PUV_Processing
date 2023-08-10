clear
figureDir = '~/GoogleDriveUCSD/MOPS_Climatology/';
MOPname = 'D0667'; % Cardiff Seaside
% MOPname = 'D0158'; % Coronado
% MOPname = 'D0045'; % IB

date1 = datetime(2000,1,1);
date2 = datetime(2019,1,1);
MOP = read_MOPline(MOPname,date1,date2);
spec1D = MOP.spec1D;
%% get basic variables
f = MOP.frequency;
fswell = [0.04 0.25];
iswell = find(f>=fswell(1) & f<=fswell(2)); % determine swell freqs
df = double(MOP.fbounds(2,:)-MOP.fbounds(1,:))'; % f bounds are not constant!!
nt = length(MOP.time); % number of observations

%% Locate peaks, get info
spec0 = [zeros(length(MOP.spec1D),1) MOP.spec1D ]; % add a zero for lowest frequency, in case there is a peak at f=0.04.
specnorm = spec0./ max(spec0,[],2); % normalize the spectrum for finding peaks

% preallocate variables for the loop
fpeaks = nan(3,nt); % <-- assumes there are never more than 3 peaks
numpks = nan(1,nt);
A(nt).loc = 0; % preallocate structure that holds peak information

for i=1:nt
    [pks,loc,w,p] = findpeaks(specnorm(i,:)); % using normalized spectrum
    
    sigpeak = p>0.5; % if a peak is not at least half the main peak, disregard.
    
    A(i).pks = pks(sigpeak); % peak height (normalized)
    A(i).loc = loc(sigpeak); % indices of peak(s) location 
    A(i).w = w(sigpeak); % peak width
    A(i).p = p(sigpeak); % peak prominence
    numpks(i) = numel(pks(sigpeak)); % number of peaks found
    
    fpeaks(1:numpks(i),i) = f(A(i).loc); % track bimodality

end
%% get bulk parameters

% peak frequency (direct from MOP model)
fpeak = MOP.fp;

% centroid frequency, spread and kurtosis
[fcentroid_swell, fspread, fkurt] = get_fstats(f(iswell),MOP.spec1D(:,iswell),df(iswell));
[fcentroid_swell, fspread, fkurt] = get_fstats(f,MOP.spec1D,df);


% total Energy
Etot = sum ( MOP.spec1D .* df',2 );

[EfluxX, EfluxY] = get_EfluxMOP(MOP);

% lowest peak
fpeaklowest = nan(size(fpeak));
for i=1:nt
    fpeaklowest(i) = f(A(i).loc(1));
end

% sqrt( Ho Lo)
Tp = 1./MOP.fp;
Lo = 9.8/(2*pi)*Tp.^2;
Ho = reverseshoal(MOP.Hs,10,MOP.fp);

sqrtHoLo = sqrt(Ho.*Lo);

% find indices of bimodal (pk2) spectra
bimodal = (numpks == 2);

monthnum = MOP.time.Month;

% months departure from August, this is beyond klugey but whatever
for i=1:nt
    if monthnum(i)>2
        mshift(i)  = abs(monthnum(i)-8);
    elseif monthnum(i)<=2
        mshift(i) = monthnum(i)+4;
    end    
end

%% make plot
close all
Z.label = 'Mean f_p (Hz)';
Z.labelInterpreter = 'tex';
Z.array = fpeak;
Z.binMethod = 'mean';
Z.cmap = cmocean('speed');
Z.scale = 'linear';
% Z.label = 'Fraction bimodal';
% Z.labelInterpreter = 'tex';
% Z.array = bimodal;
% Z.binMethod = 'fraction';


% Z.label ='$\sqrt (HoLo)$';
% Z.labelInterpreter = 'latex';
% Z.array = sqrtHoLo;
% Z.binMethod = 'mean';

Y.label = '$\sqrt{HoLo}$ (m)';
Y.labelInterpreter = 'latex';
Y.array = sqrtHoLo;
% X.bin = [245:3:295];
% X.bin = [0.25:0.25:5.25];
% X.bin = [7.75:0.25:15];
% X.bin = f;
Y.bin = [0:1:50];
Y.scale = 'linear';

% Y.label = 'E Flux';
% Y.labelInterpreter = 'tex';
% Y.array = Eflux;
% Y.bin = [1000:2000:140000];
% Y.scale = 'log';


X.label = 'Frequency Spread (Hz)';
X.labelInterpreter = 'tex';
X.array = fspread;
X.bin = [0.01:0.005:0.12]; % I made these bounds up
X.scale = 'linear';

HC.xtick = [0 1];
HC.xticklabel= {'summer','winter'};
HC.label =  'Avg. Season';
HC.array = mshift;

% HC.xtick = [0 1];
% HC.xticklabel= {'0','1'};
HC.label =  'Fraction bimodal';
HC.array = bimodal;
HC.binMethod = 'fraction';
% HC.yLimit = 1000;

figureDir = '~/GoogleDriveUCSD/MOPS_Climatology/';
figureName = [figureDir 'sqrtHoLo_fspread_bimodal_IB_2.jpeg'];

% [zdata, hdata, ZBIN]
[bimodalmat,hdatabimodal, ZBIN] = plot_jointpdfMOPS(X,Y,Z,HC,figureName,0);


% %% Pick out spectra for testing
% close all
% figure
% tol = eps(0.5);
% xind = find(abs(ZBIN.xind-0.045) < tol); % find fspread in bin = 0.045-0.05 Hz
% yind = find(abs(ZBIN.yind-34) < tol); % find sqrt(HoLo)in bin = 33-34
% zind = ZBIN.zind{xind,yind}; % get indices for spectra in this bin
% clf
% % cmap = cmocean('thermal',numel(zind));   
% cmap = parula(numel(zind));
% % for i=1:length(zind)
% % i = 6
% % hp = plot(MOP.frequency,MOP.spec1D(zind([6 11]),:),'linewidth',2) ;
% % hl = legend(datestr(MOP.time(zind([6 11]))));
% for i=1:numel(zind)
%     hp(i) = plot(MOP.frequency,MOP.spec1D(zind(i),:),'linewidth',2,'color',cmap(i,:)) ;
%     hold on
% end
% hl = legend(datestr(MOP.time(zind)));
% hl.Box = 'off';
% hl.FontSize = 14;
% colorlegend(hl,hp)
% title('$\sqrt(H_0L_0)$ = 33-34 m, fspread = 0.045-0.05 Hz', 'interpreter','latex')
% xlabel('Energy (m^2/Hz)')
% ylabel('Frequency (Hz)')
% % pause
% 
% % end   
% 
% t = MOP.time(zind([6 11]));
% % t = MOP.time(zind)
% 
% % run_makeInputFilesMOPS(t)
% % figureName = [figureDir 'specCASES_34_0p045.jpeg'];
% % print(gcf, '-djpeg', figureName,'-r300');