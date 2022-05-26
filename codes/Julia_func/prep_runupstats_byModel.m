%%
clear
filedir1 = '~/Documents/NeoStockdon/data/';
load([filedir1 'binnedMOPS.mat'],'ZBIN','fspread','fcentroid_swell','MOP','fkurt')

% xind = ZBIN.xind;
%
SigS06 = [];
Sig = [];
Eig = [];
Sinc = [];
SincS06 = [];
etaRunup = [];
etaRunupS06 = [];
R2 = [];
R2S06 = [];
sqrtHoLo = [];
fsp = [];
H = [];
fm = [];
kurt = [];
fp = [];
tIndex = [];
beta = [];

for j=25:47
    filedir2 = ['~/Documents/NeoStockdon/processed/D0045/' num2str(j,'%02.0f') '/'];

        try
            IGswash = load([filedir2 'runupstats.mat']);
            nBin = length(IGswash.Sig);
            SigS06tmp = IGswash.Sig_S06;
            sqrtHoLotmp = IGswash.sqrtHoLo;
            Sigtmp = IGswash.Sig;
            Eigtmp = IGswash.Eig;
            
            Sinctmp = IGswash.Sinc;
            SincS06tmp = IGswash.Sinc_S06;

            etaRunuptmp = IGswash.etaRunup;
            etaRunupS06tmp = IGswash.etaRunup_S06;
            
            R2tmp = IGswash.R2_stat;
            R2S06tmp = IGswash.R2_S06;

            Sigtmp = IGswash.Sig;

            tind = IGswash.tind;
            Htmp = IGswash.Ho;
            btmp = IGswash.beta; 
        end
        
        if exist('R2tmp','var')
            SigS06 = [SigS06 SigS06tmp];
            SincS06 = [SincS06 SincS06tmp];
            etaRunupS06 = [etaRunupS06 etaRunupS06tmp];
            R2S06 = [R2S06 R2S06tmp];

            sqrtHoLo = [sqrtHoLo sqrtHoLotmp];
            Sig = [Sig Sigtmp];
            Eig = [Eig Eigtmp];
            
            Sinc = [Sinc Sinctmp];
            etaRunup = [etaRunup etaRunuptmp];
            R2 = [R2 R2tmp];

            fsp = [fsp fspread(tind)];
            fm = [fm fcentroid_swell(tind)];
            H = [H Htmp];
            fp = [fp MOP.fp(tind)'];
            kurt = [kurt fkurt(tind)];
            tIndex = [tIndex tind];
            beta = [beta btmp];
            
            clear *tmp
        end
end

tdate = MOP.time(tIndex);
clear *tmp* tindBoth tindEig0 fspread i* fcentroid_swell j MOP ZBIN fkurt  filedir* tind nBin IGswash
%%

t = load('~/Documents/NeoStockdon/processed/runupstats_Bound15.mat','tIndex');
if length(tdate)>length(t.tIndex)
    [~, ~, ib] = intersect(t.tIndex,tIndex,'stable');
    Eig         = Eig        (ib);                
    H           = H          (ib);                
    R2          = R2         (ib);                
    R2S06       = R2S06      (ib);                
    Sig         = Sig        (ib);                
    SigS06      = SigS06     (ib);                
    Sinc        = Sinc       (ib);                
    SincS06     = SincS06    (ib);                
    beta        = beta       (ib);                
    etaRunup    = etaRunup   (ib);                
    etaRunupS06 = etaRunupS06(ib);                
    fm          = fm         (ib);                
    fp          = fp         (ib);                
    fsp         = fsp        (ib);                
    kurt        = kurt       (ib);                
    sqrtHoLo    = sqrtHoLo   (ib);                
    tIndex      = tIndex     (ib);                
    tdate       = tdate      (ib);    
    clear ib t
end
save '~/Documents/NeoStockdon/processed/runupstats_Bound.mat'