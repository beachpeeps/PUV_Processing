%% PROCESS PUV DATA
% need data as hourly intervals


depth_all = nanmean(PUV.P) + PUV.offsets.doffp;  % This is the hourly averaged water depth
tide_pred = detrend(depth_all); % the offset is weird

doffp = PUV.offsets.doffp;
%%
for kk = 1:1%size(PUV.time,2)
    time = PUV.time(:,kk);
    P = PUV.P(:,kk);
    U = PUV.U(:,kk);
    V = PUV.V(:,kk);
    W = PUV.W(:,kk);
    depth = depth_all(kk)
    
    fs = PUV.fs;
    dt = 1/fs;
    

    U = detrend(U);
    V = detrend(V);
    P = detrend(P);
    
    %% CHECK WHAT THIS IS
    p=15;
    N = length(U);
    
    % this will produce 2P-1 tapers, so 13
    [PSI,~]=sleptap(N,p);
    %%

    [~,Suu,Svv,Suv]=mspec(dt,U,V,PSI);
    [~,Spp,~,Spu]=mspec(dt,P,U,PSI);
    [fm,~,~,Spv]=mspec(dt,P,V,PSI);

    df = fm(2)-fm(1);
    
    %%
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
correction(ii) = cosh(k(ii)*depth) ./ cosh(k(ii)*doffp);
convert(ii) = (2*pi*fm(ii)./(g*k(ii))) .* cosh(k(ii)*doffp) ./ cosh(k(ii)*doffu);  % record conversion coefs to check them
%units = m*s^2/(m*s) = s


SSE = Spp .* (correction.^2) ;   % correct pressure for attentuation to get sea-surface elevation
SSEt = SSE(ii);

Suut = Suu(ii);
Svvt = Svv(ii);
Suvt = Suv(ii);

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
%% SPECTRAL WEIGHTED AVERAGES & STATS
% find indices of freq bands

i_ig = find(fm>0.004 & fm<0.04);
i_iglo = find(fm>=0.004 & fm<=0.025);
i_iglo = find(fm>0 & fm<=0.025);
i_ighi = find(fm>=0.025 & fm<=0.04);
i_swell = find(fm>0.04 & fm<0.2);
i_swelllo = find(fm>0.04 & fm<0.1);
i_sea = find(fm>0.1 & fm<0.25);
i_lo = find(fm>0 & fm<0.004);

Ztest =  Spp(i_swell) ./  SppU(i_swell);
ztest_all = Spp./SppU;

ztest_ss =  sum ( (Spp(i_swell) ./  SppU(i_swell) ).* SSE(i_swell) )/sum(SSE(i_swell));
%% SHEAR WAVE VELOCITY
% find shear wave velocity using R method, after Lippmann et al. 1999
q2F = Suu+Svv; %eq. 3a Noyes et al,2001. After Lippmann (1999)
r = q2F*depth/g./Spp;

q2swF = q2F-g/depth.*Spp; %shear wave velocity variance, Noyes eq. 4
q2swF(r<=1) = 0;

q2sw = sum(q2swF(i_ighi)).*df; %shear wave velocity variance, Noyes eq. 4
q2 = sum(q2F(i_ighi)).*df; %iG wave velocity variance, Noyes eq. 4



q2swlo = sum(q2swF(i_lo)).*df; %shear wave velocity variance, Noyes eq. 4
q2lo = sum(q2F(i_lo)).*df; %lo IG wave velocity variance, Noyes eq. 4
n2 = sum(Spp(i_ighi)).*df;
n2lo = sum(Spp(i_lo)).*df;
lippmann = q2*depth/g/n2;

%% COHERENCE & PHASE, etc
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

%%  WAVE DIRECTION & SPREAD
%  before rotation 0 deg is for waves headed towards positive x (usually EAST)

a1 = coPUpres ./ sqrt( Spp .* ( UUpres + VVpres ) );
b1 = coPVpres ./ sqrt( Spp .* ( UUpres + VVpres ) );
dir1 = radtodeg ( atan2(b1,a1) );
spread1 = radtodeg ( sqrt( 2 .* ( 1-sqrt(a1.^2 + b1.^2) ) ) );
spread1_ss = nanmean(spread1(i_swell));
spread1_ig = nanmean(spread1(i_ig));

a2 = (UUpres - VVpres) ./ (UUpres + VVpres);
b2 = 2 .* coUVpres ./ ( UUpres + VVpres );
dir2 = radtodeg ( atan2(b2,a2)/2 );
spread2 = (180/pi) .* sqrt(abs( 0.5 - 0.5 .* ( a2.*cos(2.*degtorad(dir2)) + b2.*sin(2.*degtorad(dir2)) )  ));
spread2_ss = nanmean(spread2(i_swell));
spread2_ig = nanmean(spread2(i_ig));

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

%% RADIATION STRESS ESTIMATES
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

Sxx = ( (1.5 + 0.5*a2) .* (Cg./C) - 0.5 ) .* SSE;
Syy = ( (1.5 - 0.5*a2) .* (Cg./C) - 0.5 ) .* SSE;
Sxy = 0.5*b2 .* (Cg./C) .* SSE;

%%% Sea Swell Bands
% These can be defined just for the swell bands
Sxx_ss = sum(Sxx(i_swell)) * df;
Sxy_ss = sum(Sxy(i_swell)) * df;


Eflux_ss = sum(Eflux(:,i_swell),2) * df;
Eflux_ig = sum(Eflux(:,i_ig),2) * df;
Eflux_ighi = sum(Eflux(:,i_ighi),2) * df;
Eflux_iglo = sum(Eflux(:,i_iglo),2) * df;
Eflux_sw = sum(Eflux(:,i_swelllo),2)*df;
Eflux_sea = sum(Eflux(:,i_sea),2)*df;

Cpu_ss = sum( cohPUpres(i_swell).*SSE(i_swell) ) / sum(SSE(i_swell));

%%
df = fm(2)-fm(1);
% significant wave height (m):

Hsig_ss = 4 * sqrt( sum( SSE(i_swell) * df ) );
Hsig_ig = 4 * sqrt( sum( SSE(i_ig) * df ) );
Hsig_iglo = 4 * sqrt( sum( SSE(i_iglo) * df ) );
Hsig_ighi = 4 * sqrt( sum( SSE(i_ighi) * df ) );
Hsig = 4 * sqrt( nansum( SSE * df ) );

% UVsig_ss = 4 * sqrt( sum( vel(i_swell) * df ) );
% UVsig_ig = 4 * sqrt( sum( vel(i_ig) * df ) );
UVsig_ss = nan; UVsig_ig = nan;

Usig_ss = 4 * sqrt( sum( Suu(i_swell) * df ) );
Usig_ig = 4 * sqrt( sum( Suu(i_ig) * df ) );
Usig_ighi = 4 * sqrt( sum( Suu(i_ighi) * df ) );


Vsig_ss = 4 * sqrt( sum( Svv(i_swell) * df ) );
Vsig_ig = 4 * sqrt( sum( Svv(i_ig) * df ) );
Vsig_ighi = 4 * sqrt( sum( Svv(i_ighi) * df ) );

%Hsig_ig = 4 * sqrt( sum( SSE(i_ig) * df ) );
%Hsig = Hsig_swell;

%HsigU_swell = 4 * sqrt( sum( SSEU(i_swell) * df ) );
%HsigU_ig = 4 * sqrt( sum( SSEU(i_ig) * df ) );
%HsigU = HsigU_swell;
% swell direction:

%dir1_swell = sum( dir1(i_swell).*SSEU(i_swell) ) / sum(SSEU(i_swell)) ;
%spread1_swell = sum( spread1(i_swell).* SSEU(i_swell) ) / sum(SSEU(i_swell)) ;
%dir2_swell = sum( dir2(i_swell).* SSEU(i_swell) ) / sum(SSEU(i_swell)) ;
%spread2_swell = sum( spread2(i_swell).* SSEU(i_swell) ) / sum(SSEU(i_swell)) ;


% dir1 is order 1, dir 2 is order 2
dirt  = [ dir1(ii);  dir2(ii)]  ;
spreadt =  spread2(ii)  ;

a1t = a1(ii);
a2t = a2(ii);
b1t = b1(ii);
b2t = b2(ii);

a1_ss = sum(  a1(i_swell).*SSE(i_swell) ) / sum(SSE(i_swell)) ;
b1_ss = sum(  b1(i_swell).*SSE(i_swell) ) / sum(SSE(i_swell)) ;
a2_ss = sum(  a2(i_swell).*SSE(i_swell) ) / sum(SSE(i_swell)) ;
b2_ss = sum(  b2(i_swell).*SSE(i_swell) ) / sum(SSE(i_swell)) ;


dir_ss = rad2deg ( atan2(b2_ss,a2_ss)/2 );     %based off u and v alone (not dependent on p!)
spread_ss = (180/pi) .* sqrt(abs( 0.5 - 0.5 .* ( a2_ss.*cos(2.*deg2rad(dir_ss)) + b2_ss.*sin(2.*deg2rad(dir_ss)) )  ));

dir1_ss = rad2deg(atan2(b1_ss,a1_ss));
% spread1_ss = (180/pi).*sqrt(2*(1-sqrt(a1_ss^2 + b1_ss^2)));

% centroid frequency
fcentroid_swell = sum ( fm(i_swell).* SSE(i_swell) ) / sum ( SSE(i_swell) ) ;
Tm = 1/fcentroid_swell;

% frequency spread
fspread = sum ( (fm(i_swell)-fcentroid_swell).^2.* SSE(i_swell) ) / sum ( SSE(i_swell) ) ;



% peak frequency
i_peak = find(SSE == max(SSE));
fpeak = fm(i_peak);

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
