function get_data_sectormoments(datahour,savename,casenumber,inorm)
% clear
% datahour = datenum(2013,10,21,13,0,0);
% % datahour = datenum(2013,10,1,22,0,0);
% 
% savename = '../mat/moments_data_IN_22';
% casen = 22;

casen = casenumber; 
%%
load('~/Agate/mat/sensorlocs.mat', 'X','indP');
xloc = X(indP);
Hz = 2;
Pss = nan(7167,18);
currdir = pwd;

cd ~/Agate/mfiles
for i=1:18
    try
        [H,Pss(:,i),time(:,i),Zstart,Zsurvey,pressure] = pcorrect_agate(i,datahour,Hz);
    catch
    end
end
cd(currdir)
%%
%%
cd ~/Agate/mfiles/PUVprocessing/
[xloc_PUV, PssPUV, numvec] = get_hourlyPatPUV(datahour);
cd(currdir)
% PssPUV(:,3) = [];
% xloc_PUV(3)=[];
%%
X = [xloc(:); xloc_PUV(:)];
[B,I] = sort(X);
P = [Pss PssPUV];
P = P(:,I);
Pss = P;
xloc = B;

%%
load ~/Agate/mat/swashbase.mat
tvec = swashbase.tvec;
X = swashbase.X;
indT = knnsearch(tvec(:),datahour);

swashloc = X(swashbase.swashind2(indT));

Ptonan = find(xloc<swashloc);
Pss(:,Ptonan) = nan(7167,length(Ptonan));

%% nan the duplicate sensors
Pss(:,[2 7 12 14]) = nan(7167,4);
xarray = xloc;
%%
%     xarray(isnan(Pss(1,:))) = [];
%     Pss(:,isnan(Pss(1,:))) = [];
%%
T = 7100/2;
sig = Pss(1:7100,:);
    flim = [0.004 0.025 0.04 0.25];
flim = [0.004 0.04 0.25];


[m2,m3] = sectormoments( T, sig , flim,inorm);
%     [m2_1,m3_1] = sectormoments( T, sig , flim,1);

%%
biphase.self = atan2(imag(m3.self),real(m3.self));
biphase.sum = atan2(imag(m3.sum),real(m3.sum));
biphase.dif = atan2(imag(m3.dif),real(m3.dif));
%%
biphase_self(casen,:,:)  = biphase.self;
biphase_sum(casen,:,:)  = biphase.sum;
biphase_dif(casen,:,:)  = biphase.dif;
%
skew_self(casen,:,:) = real(m3.self);
skew_sum(casen,:,:) = real(m3.sum);
skew_dif(casen,:,:) = real(m3.dif);
skew_tot(casen,:,:) = real(m3.tot);
% norm3(casen,:,:) = m3.norm3;

asym_self(casen,:,:) = imag(m3.self);
asym_sum(casen,:,:) = imag(m3.sum);
asym_dif(casen,:,:) = imag(m3.dif);
asym_tot(casen,:,:) = imag(m3.tot);

bicoh_self(casen,:,:) = sqrt(dot(m3.self,m3.self,3));
bicoh_sum(casen,:,:) = sqrt(dot(m3.sum,m3.sum,3));
bicoh_dif(casen,:,:) = sqrt(dot(m3.dif,m3.dif,3));


vartot(casen,:,:) = m2;


%%

clearvars -except t skew* biphase* bicoh* asym* savename xarray vartot s2 norm3 xloc_PUV xloc
clear biphase
save(savename)
disp(['Data moments written to' savename])
