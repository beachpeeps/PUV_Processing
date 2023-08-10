function [A] = vector_wave_stats_mspec(U,V,P,tt, doffp, doffu)
% addpath jlab


depth = nanmean(P) + doffp;  % This is the hourly averaged water depth

%%
%%%%%%%%%%%%%%%%%%%%%%

U = detrend(U);
V = detrend(V);
P = detrend(P);

% dsample = 1;
% P = downsample(P,dsample);
% U = downsample(U,dsample);
% V = downsample(V,dsample);

N = length(U);

% U = U(1:N)'; V = V(1:N)'; P = P(1:N)';
p = 15;
fs = 2;  % Hz
dt = 1/fs;
% this will produce 2P-1 tapers, so 13
[PSI,~]=sleptap(N,p);

if isrow(U)
    U = U';
end
if isrow(V)
    V = V';
end
if isrow(P)
    P = P';
end

[~,Suu,Svv,Suv]=mspec(dt,U,V,PSI);
[~,Spp,~,Spu]=mspec(dt,P,U,PSI);
[fm,~,~,Spv]=mspec(dt,P,V,PSI);

id_low = find(min(abs(fm - 0.004))== abs(fm-0.004));
Spp(1:id_low) = NaN; Suu(1:id_low) = NaN; Suv(1:id_low) = NaN;
Spu(1:id_low) = NaN; Spv(1:id_low) = NaN; Suv(1:id_low) = NaN;

df = fm(2)-fm(1);
time = tt;
%% DONE Spectra 
% DEPTH CORRECTION -------------------------------------------------------------
% find correction/conversion rcoefs at each f-band
% to calc sea surface elevation and convert velocities to pressure units



g = 9.81;    % gravity, m/s^2
omega = 2*pi*fm;   % This is radian frequency

k = get_wavenumber( omega , depth);


correction = zeros(size(fm));
convert = zeros(size(fm));

fcutoff = 0.25;

ii = find(fm <= fcutoff);
i_ig = find(fm>0.004 & fm<0.04);
i_swell = find(fm>0.04 & fm<0.2);


correction(ii) = cosh(k(ii)*depth) ./ cosh(k(ii)*doffp);
%limit the correction factor K_p = correction^2 to 10 (order of magnitude)
%this will allow for adjustment in depth and stop overcorrection
%attempted to follow Jones and Monismith 2007 at 12* the noise floor 
%but their sensors are too shallow. 
cutoffFreq = fm(find(correction.^2 >= 10,1));
correction(correction.^2>10) = sqrt(10);
convert(ii) = (2*pi*fm(ii)./(g*k(ii))) .* cosh(k(ii)*doffp) ./ cosh(k(ii)*doffu);  % record conversion coefs to check them
%units = m*s^2/(m*s) = s


SSE = Spp .* (correction.^2) ;   % correct pressure for attentuation to get sea-surface elevation
SSEt = SSE(ii);
Sppt = Spp(ii);

Suut = Suu(ii);
Svvt = Svv(ii);
Suvt = Suv(ii);
Sput = Spu(ii);
Spvt = Spv(ii);

fmt = fm(ii);

UUpres = Suu .* (convert.^2) ;   % convert velocities into "equivalent pressure" (to compare with PP)
VVpres = Svv .* (convert.^2) ;
UVpres = Suv .* (convert.^2) ;
%units = m^2/s^2 /Hz* s^2 = m^2*s

PUpres = Spu .* convert ;   % convert x-spectra as well (but conversion is not squared here)
PVpres = Spv .* convert ;
% units = m^2/s*Hz*s = m^2*s

% should check that:   PP^2 = UUpres^2 + VVpres^2
% the 'pres' velocities are now the right units to compare with P
SppU = UUpres+VVpres; % equivalent to p^2 = u^2+v^2
SSEU = SppU .* (correction.^2); %sea surface elevation determined from velocity
SSEUU = UUpres.*(correction.^2);%sea surface elevation determined from U velocity
SSEUP = PUpres.*correction; %sea surface elevation determined from PU cospectra

Spec = struct('fm', fmt,'SSE',SSEt, 'Spp', Sppt,'Suu', Suut, 'Svv', Svvt, 'Spu', Sput, 'Spv', Spvt, 'Suv', Suvt, ...
            'UUpres', UUpres, 'VVpres', VVpres, 'UVpres', UVpres, 'PUpres', PUpres, 'PVpres', PVpres, ...
            'SppU', SppU, 'SSEU', SSEU, 'SSEUU',SSEUU, 'SSEUP', SSEUP,'cutoffFreq',cutoffFreq);
        
%% DONE Z TEST
ztest_all = Spp./SppU;
ztest_ss =  Spp(i_swell) ./  SppU(i_swell);
ztest_ig =  Spp(i_ig) ./  SppU(i_ig);

ztest_ss_sum =  sum((Spp(i_swell)./ SppU(i_swell)).* SSE(i_swell))/sum(SSE(i_swell));
ztest_ig_sum =  sum((Spp(i_ig)./ SppU(i_ig)).* SSE(i_ig))/sum(SSE(i_ig));

ztest = struct('ztest_all', ztest_all, 'ztest_ss', ztest_ss, 'ztest_ig', ztest_ig,...
            'ztest_ss_sum', ztest_ss_sum, 'ztest_ig_sum', ztest_ig_sum);
        
%% SHEAR WAVE VELOCITY
% % find shear wave velocity using R method, after Lippmann et al. 1999
% q2F = Suu+Svv; %eq. 3a Noyes et al,2001. After Lippmann (1999)
% r = q2F*depth/g./Spp;
% 
% q2swF = q2F-g/depth.*Spp; %shear wave velocity variance, Noyes eq. 4
% q2swF(r<=1) = 0;
% 
% q2sw = sum(q2swF(i_ighi)).*df; %shear wave velocity variance, Noyes eq. 4
% q2 = sum(q2F(i_ighi)).*df; %iG wave velocity variance, Noyes eq. 4
% 
% 
% 
% q2swlo = sum(q2swF(i_lo)).*df; %shear wave velocity variance, Noyes eq. 4
% q2lo = sum(q2F(i_lo)).*df; %lo IG wave velocity variance, Noyes eq. 4
% n2 = sum(Spp(i_ighi)).*df;
% n2lo = sum(Spp(i_lo)).*df;
% lippmann = q2*depth/g/n2;

%% DONE COHERENCE & PHASE, etc
% Definitions from Elgar, 1985, JFM ... gives same result as Wunsch's def.
% using pressure converted or regular velocity should give the same coh & phase
% choose to use pressure converted in order to calc dirs later
% Cospectrum & Quadrature:

coPUpres = real(PUpres);   quPUpres = imag(PUpres);
coPVpres = real(PVpres);   quPVpres = imag(PVpres);
coUVpres = real(UVpres);   quUVpres = imag(UVpres);

%%% Coherence & Phase at each freq-band
% *** note that it's important to calc this AFTER all merging and ensemble avg.

warning off MATLAB:divideByZero
cohPUpres = sqrt( (coPUpres.^2 + quPUpres.^2) ./ (Spp.* UUpres) );
phPUpres  = (180/pi) .* atan2( quPUpres , coPUpres );
cohPVpres = sqrt((coPVpres.^2 + quPVpres.^2)./ (Spp.* VVpres));
phPVpres  = (180/pi) .* atan2( quPVpres , coPVpres );
cohUVpres = sqrt((coUVpres.^2 + quUVpres.^2)./(UUpres .* VVpres));
phUVpres  = (180/pi) .* atan2( quUVpres , coUVpres );

coh = struct('cohPUpres', cohPUpres, 'phPUpres', phPUpres, 'cohPVpres', cohPVpres, 'phPVpres', phPVpres, 'cohUVpres', cohUVpres, 'phUVpres', phUVpres);

%% DONE FOURIER COEFFICIENTS
%  before rotation 0 deg is for waves headed towards positive x (usually EAST)
% Herbers, Elgar and Guza 1999 - Directional spreading of waves in the nearshore

a1 = coPUpres ./ sqrt( Spp .* ( UUpres + VVpres ) );
b1 = coPVpres ./ sqrt( Spp .* ( UUpres + VVpres ) );
a2 = (UUpres - VVpres) ./ (UUpres + VVpres);
b2 = 2 .* coUVpres ./ ( UUpres + VVpres );

a1t = a1(ii);
a2t = a2(ii);
b1t = b1(ii);
b2t = b2(ii);

a1_ss = sum( a1t(i_swell).*SSE(i_swell) ) / sum(SSE(i_swell)) ;
b1_ss = sum( b1t(i_swell).*SSE(i_swell) ) / sum(SSE(i_swell)) ;
a2_ss = sum( a2t(i_swell).*SSE(i_swell) ) / sum(SSE(i_swell)) ;
b2_ss = sum( b2t(i_swell).*SSE(i_swell) ) / sum(SSE(i_swell)) ;

FC = struct('a1', a1t, 'b1', b1t, 'a2', a2t, 'b2', b2t, 'a1_ss', a1_ss, 'b1_ss', b1_ss, 'a2_ss', a2_ss, 'b2_ss', b2_ss);

%% DONE WAVE DIRECTION & SPREAD
dir1 = radtodeg ( atan2(b1t,a1t) );
spread1 = radtodeg ( sqrt( 2 .* ( 1-sqrt(a1t.^2 + b1t.^2) ) ) );
spread1_ss = nanmean(spread1(i_swell));
spread1_ig = nanmean(spread1(i_ig));

dir2 = radtodeg ( atan2(b2t,a2t)/2 );
spread2 = (180/pi) .* sqrt(abs( 0.5 - 0.5 .* ( a2t.*cos(2.*degtorad(dir2)) + b2t.*sin(2.*degtorad(dir2)) )  ));
spread2_ss = nanmean(spread2(i_swell));
spread2_ig = nanmean(spread2(i_ig));

dir1_ss_sum = rad2deg(atan2(b1_ss,a1_ss));
spread1_ss_sum = (180/pi).*sqrt(2*(1-sqrt(a1_ss^2 + b1_ss^2)));

dir2_ss_sum = rad2deg ( atan2(b2_ss,a2_ss)/2 );     %based off u and v alone (not dependent on p!)
spread2_ss_sum = (180/pi) .* sqrt(abs( 0.5 - 0.5 .* ( a2_ss.*cos(2.*deg2rad(dir2_ss_sum)) + b2_ss.*sin(2.*deg2rad(dir2_ss_sum)) )  ));

dir = struct('dir1', dir1, 'spread1', spread1, 'spread1_ss', spread1_ss, 'spread1_ig', spread1_ig, 'dir1_ss_sum', dir1_ss_sum, 'spread1_ss_sum', spread1_ss_sum,...
             'dir2', dir2, 'spread2', spread2, 'spread2_ss', spread2_ss, 'spread2_ig', spread2_ig, 'dir2_ss_sum', dir2_ss_sum, 'spread2_ss_sum', spread2_ss_sum);
    
%% ESTIMATE ENERGY FLUX
% spectral estimator from T.H.C. Herbers, using group velocity:

C = omega./k;
Cg = get_cg(k,depth);   %0.5 * omega.* k .* ( 1 + (2*k*depth)./sinh(2*k*depth) );
rho = 1028;    % density (Kg /m^-3)
const = rho*g * Cg.* ((cosh(k*depth)).^2) ./ ((cosh(k*doffp)).^2);
%units = kg/m^3* m/s^2*m/s = kg/ms^3
om = omega;
CO = const;

%%% Energy Flux (by freq) in cartesian coordinates
% From Sheremet et al (2002),
% Eflux (pos) = 1/4*Co_pp + 1/4*(h/g)*Co_uu + 1/2*sqrt(h/g)*Co_pu
% This assumes shallow water and cross-shore propagation, so we alter this
% to have frequency dependent group speed (in the constant), and UUpres and
% VVpres etc are already in terms of sea surface elevation (units: meters^2/Hz),
% so no need for the (h/g) or sqrt(h/g) terms. Thus, it changes into:

posX=const.*(0.25.* ( abs(Spp) + (UUpres) ) + 0.5*real(PUpres));
negX=const.*(0.25.* ( abs(Spp) + (UUpres) ) - 0.5*real(PUpres));
posY=const.*(0.25.* ( abs(Spp) + (VVpres) ) + 0.5*real(PVpres));
negY=const.*(0.25.* ( abs(Spp) + (VVpres) ) - 0.5*real(PVpres));
%units = kg/ms^3 * m^2*s = kg*m/s^2
% when this gets integrated over the ss or ig frequencies it becomes
% kg*m/s^3: take out rho (kg/m^3) and you get m^4/s^3

posX2= rho*g * Cg .* a1 .* SSE;
posY2 = rho*g * Cg .* b1 .* SSE;
% units = kg/m^3 * m/s^2 * m/s * m^2 *s = kg*m/s^2

% *** POSSIBLE PROBLEM WITH NEGATIVE ENERGY FLUXES ***
% use 2 lines below to set negative energy fluxes to be very small positives:
% posX( posX<10^-7 ) = 10^-7; negX( negX<10^-7 ) = 10^-7;
% posY( posY<10^-7 ) = 10^-7; negY( negY<10^-7 ) = 10^-7;
%create energy flux variable (for output)

Eflux = [ posX2 posX negX posY negY]';
Epos = Eflux(2,:)*df; Eneg = Eflux(3,:)*df;

Eflux_dir = struct('posX', posX(ii), 'posY', posY(ii), 'negX', negX(ii), 'negY', negY(ii), 'Cg', Cg(ii), 'depth', depth);

%% DONE RADIATION STRESS ESTIMATES
% Here following eqns 1a-1c in Herbers and Guza (1989), after Longuet
% Higgins
%
% $$S_{xx}(f) = \rho g \int_0^{2\pi}d\theta E(f,\theta)[n(f) -1/2+n(f)\cos^2\theta]$$
%
% $$S_{xy}(f) = S_{yx}(f) = \rho g \int_0^{2\pi}d\theta E(f,\theta)n(f)\cos\theta\sin\theta$$
%
% $$S_{yy}(f) = \rho g \int_0^{2\pi}d\theta E(f,\theta)[n(f) - 1/2+n(f)\sin^2\theta]$$
%
% where n is the ratio of group speed to phase speed

c = Cg./C; ct = c(ii);
Sxx = ( (1.5 + 0.5*a2t) .* ct - 0.5 ) .* SSEt;
Syy = ( (1.5 - 0.5*a2t) .* ct - 0.5 ) .* SSEt;
Sxy = 0.5*b2t .* ct .* SSEt;

%%% Sea Swell Bands
% These can be defined just for the swell bands
Sxx_ss = sum(Sxx(i_swell)) * df;
Syy_ss = sum(Syy(i_swell)) * df;
Sxy_ss = sum(Sxy(i_swell)) * df;


Eflux_ss = sum(Eflux(:,i_swell),2) * df;
Eflux_ig = sum(Eflux(:,i_ig),2) * df;

Cpu_ss = sum( cohPUpres(i_swell).*SSE(i_swell) ) / sum(SSE(i_swell));


RS = struct('Sxx',Sxx,'Syy',Syy,'Sxy',Sxy,'Sxx_ss',Sxx_ss,'Syy_ss',Syy_ss,'Sxy_ss',Sxy_ss,...
            'Eflux_ss',Eflux_ss,'Eflux_ig',Eflux_ig,'Cpu_ss', Cpu_ss);
            
%% DONE Hsig
df = fm(2)-fm(1);
% significant wave height (m):

Hs_ss = 4 * sqrt( sum( SSE(i_swell) * df ) );
Hs_ig = 4 * sqrt( sum( SSE(i_ig) * df ) );
Hs = 4 * sqrt( nansum( SSE * df ) );

Usig_ss = 4 * sqrt( sum( Suu(i_swell) * df ) );
Usig_ig = 4 * sqrt( sum( Suu(i_ig) * df ) );

Vsig_ss = 4 * sqrt( sum( Svv(i_swell) * df ) );
Vsig_ig = 4 * sqrt( sum( Svv(i_ig) * df ) );


Hsig = struct('Hs_ss', Hs_ss, 'Hs_ig', Hs_ig, 'Hs', Hs,...
    'Usig_ss', Usig_ss, 'Usig_ig', Usig_ig, 'Vsig_ss', Vsig_ss, 'Vsig_ig', Vsig_ig);

%%
%HsigU_swell = 4 * sqrt( sum( SSEU(i_swell) * df ) );
%HsigU_ig = 4 * sqrt( sum( SSEU(i_ig) * df ) );
%HsigU = HsigU_swell;
% swell direction:

%dir1_swell = sum( dir1(i_swell).*SSEU(i_swell) ) / sum(SSEU(i_swell)) ;
%spread1_swell = sum( spread1(i_swell).* SSEU(i_swell) ) / sum(SSEU(i_swell)) ;
%dir2_swell = sum( dir2(i_swell).* SSEU(i_swell) ) / sum(SSEU(i_swell)) ;
%spread2_swell = sum( spread2(i_swell).* SSEU(i_swell) ) / sum(SSEU(i_swell)) ;

%% DONE centroid frequency
fcentroid_swell = sum ( fm(i_swell).* SSE(i_swell) ) / sum ( SSE(i_swell) ) ;
Tm = 1/fcentroid_swell;

% frequency spread
fspread = sum ( (fm(i_swell)-fcentroid_swell).^2.* SSE(i_swell) ) / sum ( SSE(i_swell) ) ;



% peak frequency
i_peak = find(SSE == max(SSE));
fpeak = fm(i_peak);

ids = struct('i_swell', i_swell, 'i_ig', i_ig, 'k', k, 'omega', omega, 'fcentroid_swell', fcentroid_swell, 'fpeak',fpeak,'Tm',Tm, 'fspread',fspread, 'doffp', doffp, 'doffu', doffu);
%% ----------------------------------------------------------------------------
A = structurize;
    function A = structurize
        A = struct('time', time, 'Spec',Spec, 'ztest', ztest,'Hsig', Hsig,'FC', FC, 'coh', coh, 'dir', dir,'RS', RS,'Eflux', Eflux_dir, 'ids', ids);
    end



    function k = get_wavenumber(omega,h)
        % Falk Feddersen (c) 2001
        %
        % function that takes the radian wave frequency and
        % a vector of depths and returns the wavenumber at that
        % depth by solving the dispersion relationship
        %
        % function k = get_wavenum(omega,h)
        % returns the wavenumber of the gravity wave
        % dispersion relation, by using newtons method
        
        % the initial guess will be the shallow water wavenumber
        
        g = 9.81;
        
        k = omega./sqrt(g*h);
        
        
        f = g*k.*tanh(k.*h) - omega.^2;
        
        while max(abs(f))>1e-10
            dfdk = g*k.*h.*(sech(k.*h)).^2 + g*tanh(k.*h);
            k = k - f./dfdk;
            f = g*k.*tanh(k.*h) - omega.^2;
        end
    end



    function cg = get_cg(k,h)
        % Falk Feddersen (c) 2001
        %
        % function that takes the vector wavenumbers and depths
        % and returns the group velocity
        %
        % function cg = funwaveC_get_cg(k,h)
        g = 9.81;
        
        % check to make sure k & h are the same size
        
        [nk,mk] = size(k);
        [nh,mh] = size(h);
        
        if (nh==1)
            if (mh==1)
                cg = 0.5*(g*tanh(k*h)+g*(k.*h).*(sech(k.*h).^2))./sqrt(g*k.*tanh(k.*h));
                return;
            end;
        end;
        
        
        if (nk ~= nh)
            disp('** Error k & h have wrong sizes');
            return;
        end;
        
        if (mk ~= mh)
            disp('** Error k & h have wrong sizes');
            return;
        end;
        
        
        
        
        cg = 0.5*(g*tanh(k.*h)+g*(k.*h).*(sech(k.*h).^2))./sqrt(g*k.*tanh(k.*h));
    end

                        
end
