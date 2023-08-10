function [f , F , dF_dx , S , E] = energyBalance( t , x , u , w , p , eta , d , tstart , tend , varargin )
%
options.type = 'NH';
options.df   = 0.01;
options.fs   = 1/mean(diff(t));
options.g    = 9.81;
options.DepthAveraged = false;
options.ustar = [];
options.omega = [];
options.includeVerticalGradientTerms = false;
options.drylim = 0.1; %NOTE: in the swash zone the model sets fric velocity to 1 if not computed; this can give very weird results.
                      %      here we set all friction velocities to 0 when
                      %      the waterdepth is below drylim

disp( '____________________________________________________________________')
disp( '           Estimating Balance        ')
disp( '')
fprintf( 'Balance type: %s \n', options.type );
if options.DepthAveraged 
fprintf( 'Depth Averaged: True \n' );
else
fprintf( 'Depth Averaged: False \n' );
end
fprintf( 'number of sigma layers: %i \n', size( u , 3 ) )
fprintf( 'number of grid points: %i \n' , size(u,1));
fprintf( 'number of time steps:  %i \n' , size(u,2));
%
% If either p or w are empty arrays, do a NLSWE balance only
if isempty(w) || isempty(p)
    options.type = 'NLSWE';
end
options = parseOptions( options , varargin );
%
%If no ustar is provided, just create an array of appropriate size filled
%with 0's.
%
if isempty( options.omega )
    omega = 0*u;
    options.includeVerticalGradientTerms = false;
else
    omega = options.omega;
    options.includeVerticalGradientTerms = true;
end

if isempty( options.ustar )
    ustar = 0*eta;
else
    ustar = options.ustar;
end

%
% Do some pre-processing to make calculations easier
%
disp( 'Processing Variables')
[t,u,U,w,p,omega,eta, d,D,d_dx,d_dz,ub,tau,Qd ] = procVars( t,x,u,w,p,omega,eta,d,ustar,tstart,tend,options);

% 
% Add terms that pop up in both the non-hydrostatic and hydrostatic balance
% (note that hydrostatic terms may be depth dependend due to shear)
disp( 'Hydrostatic Terms')
[f , F ,  S , E] = NLSWEBalance( x ,u , eta  , d , D ,d_dx, ub,tau ,Qd,options);

%
% Add nonhydrostatic terms if available
%
if strcmpi(options.type,'NH') && ~options.DepthAveraged
   %
   disp( 'Non-Hydrostatic Terms')
   [f , F ,  S , E]  = nhterms(  f , F ,  S , E, x, u, w, p,omega,eta  , d , D , d_dx,d_dz, options);
   %
end
disp( '_____________________________________')
disp( '')

%Calculate flux gradients of all terms in the flux structure
names = fieldnames(F);
for jj=1:length(names)
    %
    name = names{jj};
    dF_dx.(name) = derivative( F.(name) , diff(x) );
    %
end
%
end

function [f , F ,  S , E]  = nhterms(  f , F ,  S , E, x, u, w, p,omega, eta  , d , D , d_dx,d_dz, options)
%ShortHands
g  = options.g;
df = options.df;
fs = options.fs;

d = ones( length(f),1 ) *  d;

%W term in energy

[~,Dw_w]     = cospec( D.*w, w , df , fs);
E   = E + Dw_w/2;

%Linear flux
[~,p_u]     = cospec( p, u , df , fs);
F.lin_nh = p_u .* d;



%Nonlinear flux
[~,p_uEta] = cospec( p  , eta.*u , df , fs);
[~,w_uwD] =  cospec( w  , u.*w.*D , df , fs);

%nonlinear flux
F.nl_nh= p_uEta + w_uwD / 2;
F.tot  = F.tot + F.lin_nh + F.nl_nh;

%Sources
[~,Duw_dwdx] = cospec( D.*u.*w  , d_dx.w , df , fs);
[~,Dw_udwdx] = cospec( D.*w  , u.*d_dx.w , df , fs);

[~,uEta_dpdx] = cospec( eta.*u , d_dx.p , df , fs);
[~,u_etadpdx] = cospec( u , eta .*d_dx.p  , df , fs);

%More Vertical Sources
[~,tmp]    = cospec( w    , d_dz.p , df , fs);
S.nhv   = tmp/2;
[~,tmp]  = cospec( w.*D , d_dz.p./D , df , fs);
S.nhv   = S.nhv - tmp/2;


[~,tmp]  = cospec( d_dx.s .* d_dz.p , u , df , fs);
S.nhv   = S.nhv + tmp/4;
[~,tmp]  = cospec( d_dx.s .* d_dz.p./D , u.*D , df , fs);
S.nhv   = S.nhv + tmp/4;
[~,tmp]  = cospec( d_dz.p , d_dx.s.*u , df , fs);
S.nhv   = S.nhv - tmp/2;

[~,tmp]    = cospec( omega.*u  , d_dz.u , df , fs);
S.nhv   = S.nhv + tmp/2;
[~,tmp]  = cospec( u.*D , omega.*d_dz.u./D , df , fs);
S.nhv   = S.nhv - tmp/2;

[~,tmp]    = cospec( omega.*w  , d_dz.w , df , fs);
S.nhv   = S.nhv + tmp/2;
[~,tmp]  = cospec( w.*D , omega.*d_dz.w./D , df , fs);
S.nhv   = S.nhv - tmp/2;

%source terms
S.nhp   = (uEta_dpdx - u_etadpdx)/2;
S.nhw   = Duw_dwdx/2 - Dw_udwdx/2;
S.tot   = S.tot+S.nhp+S.nhw + S.nhv;

end


function [f , F ,  S , E]  = NLSWEBalance(  x, u, eta  , d , D , d_dx,ub,tau ,Qd, options)
%ShortHands
g  = options.g;
df = options.df;
fs = options.fs;

%
% cospectra - used in energy
%
if ndims(eta)==3
    eta0 = squeeze(eta(:,1,:));
    D0 = squeeze(D(:,1,:));
else
    eta0 = eta;
    D0   = D;
end



if strcmpi( options.type , 'henderson')
    %
    % Henderson Balance
    %
    disp('- Henderson Balance')
    Q = u.*D;
    Sxx = u.*Q + g*eta.^2/2;
    M = u.*eta;
    
    [f,eta_eta]  = cospec( eta0 , eta0 , df , fs);
    [~,DQ_Q]     = cospec( Q, Q , df , fs);
    
    %
    % cospectra - used in fluxes
    %
    [~,u_eta]  = cospec( u  , eta , df , fs);
    [~,QSxx] = cospec( Qd  , Sxx , df , fs);
    [~,Meta] = cospec( M   , eta , df , fs);
    
    %           - used in source term
    [~,dQddxSxx ] = cospec( d_dx.Qd  , Sxx , df , fs);
    d = ones( length(f),1 ) *  d;
    
    %Energy
    E     = g * eta_eta/2 + DQ_Q/2./d;
    
    %Linear flux
    F.lin = g * d .* u_eta;
    
    %nonlinear flux
    F.nl   = QSxx + g*Meta;
    F.tot  = F.lin + F.nl;
    F.sw   = F.tot;
    
    
    %source terms
    [~,Qdtau] = cospec( Qd , tau  , df , fs);
    
    S.sw   = dQddxSxx;
    S.stress = Qdtau;    
    S.tot  = S.sw;
    
else
    disp('- New Balance')
    [f,eta_eta]  = cospec( eta0 , eta0 , df , fs);
    [~,Du_u]     = cospec( D.*u, u , df , fs);
    
    %
    % cospectra - used in fluxes
    %
    [~,u_eta]  = cospec( u  , eta , df , fs);
    [~,eta_uEta] = cospec( eta  , eta.*u , df , fs);
    [~,u_uuD] =  cospec( u  , u.*u.*D , df , fs);
    
    %           - used in source term
    [~,Duu_dudx] = cospec( D.*u.^2  , d_dx.u , df , fs);
    [~,Du_ududx] = cospec( D.*u  , u.*d_dx.u , df , fs);
    [~,uEta_detadx] = cospec( eta.*u , d_dx.eta , df , fs);
    [~,u_etadetadx] = cospec( u , eta .*d_dx.eta  , df , fs);
    [~,ubtau] = cospec( ub , tau  , df , fs);        
    [~,ubtauD0] = cospec( ub .* D0 , tau./D0  , df , fs);
    
    d = ones( length(f),1 ) *  d;
    
    %Energy
    E     = g * eta_eta/2 + Du_u/2;
    
    %Linear flux
    F.lin = g * d .* u_eta;
    
    %nonlinear flux
    F.nl   = g * eta_uEta + u_uuD / 2;
    F.tot  = F.lin + F.nl;
    
    %source terms
    S.sw   = Duu_dudx/2 - Du_ududx/2 + g*(uEta_detadx - u_etadetadx)/2;
    S.stress = +ubtau/2+ubtauD0/2;
    %
end







F.tot  = F.lin + F.nl;
F.sw   = F.tot;
S.tot  = S.sw+S.stress;
%

end


function [t,u,U,w,p,omega,eta,d,D, d_dx,d_dz,ub, tau,Qd ] = procVars( t,x,u,w,p,omega,eta,d,ustar, tstart,tend,options)
eta = squeeze(eta);
u   = squeeze(u);
x   = squeeze(x);

nx = size(eta,1);
nt = size(eta,2);
% Filter for time domain
if isnan(tstart)
    its=1;
else
    [~,its] = min(abs(t-tstart));
end

%Calculate derivatives
if isnan(tend)
    ite=nt;
else
    [~,ite] = min(abs(t-tend));
end

msk = its : 1 : ite;


if ndims( u ) == 3 && ~options.DepthAveraged 
    ns = size( u , 3);
    U  = sum( u , 3 ) / ns;    
    U  = U';
    U  = U(msk,:,:);
    u   = permute(u , [2,3,1]);
    u   = u(msk,:,:);
    eta = eta';
    D   = eta + ones(nt,1) * d;
    
    %note; as we multiply the free-surface with velocities it is convinient
    %to make eta(nt,nx) and array of the same size as u; i.e.
    %eta(nt,nsigma,nx). The same holds for D
    eta = repmat(eta , [1,1,ns] );
    eta = permute(eta , [1,3,2] );
    eta = eta(msk,:,:);
    
    depmsk = D < options.drylim;
    D = repmat(D , [1,1,ns] );
    D = permute(D , [1,3,2] );
    D = D(msk,:,:);
    depmsk = depmsk( msk , :);
    
    if ~isempty(p)
        %
        p = permute( p , [2,3,1]);
        w = permute( w , [2,3,1]);   
        p = p(msk,:,:);
        w = w(msk,:,:);
        omega = permute( omega , [2,3,1]);
   
        omega = omega(msk,:,:);
        wm     = zeros( size(p) );
        pm     = zeros( size(p) );
        omegam = zeros( size(p) );
        for ii = 1 : size(p,2)
            %
            wm(:,ii,:) = w(:,ii,:)/2+w(:,ii+1,:)/2;
            omegam(:,ii,:) = omega(:,ii,:)/2+omega(:,ii+1,:)/2;
            if ii == 1
                pm(:,ii,:) = p(:,ii,:)/2;
            else
                pm(:,ii,:) = p(:,ii,:)/2+p(:,ii-1,:)/2;
            end
            %
        end
        
        if options.includeVerticalGradientTerms
            %
            d_dz.p     = zeros( size(p) );
            d_dz.w     = zeros( size(p) );
            d_dz.u     = zeros( size(p) );
            d_dz.omega = zeros( size(p) );
            ns = size(p,2);
            ds = 1/ns;
            for ii = 1 : ns
                %
                d_dz.w(:,ii,:) = (w(:,ii,:)/2-w(:,ii+1,:)) / ds;
                d_dz.omega(:,ii,:) = (omega(:,ii,:)/2-omega(:,ii+1,:)) / ds;
                if ii == 1
                    d_dz.p(:,ii,:) = -p(:,ii,:)/ds;
                    d_dz.u(:,ii,:) = (u(:,ii,:)-u(:,ii+1,:) )/ds;
                else
                    d_dz.p(:,ii,:) = (-p(:,ii,:)/2+p(:,ii-1,:) )/ds ;
                    
                    if ii == ns
                       d_dz.u(:,ii,:) = (u(:,ii-1,:)-u(:,ii,:) )/ds; 
                    else
                       d_dz.u(:,ii,:) = (u(:,ii-1,:)-u(:,ii+1,:) )/ds/2.;
                    end
                end
                %
            end
            %
        end
        
        %
        w=wm;
        p=pm;
        omega=omegam;
        %w = wm(msk,:,:);
        %p = pm(msk,:,:);
        d_dx.w     = derivative( w     , diff(x) );
        d_dx.p     = derivative( p     , diff(x) );
        d_dx.omega = derivative( omega , diff(x) );
    end
    ub  = squeeze(u(:,end,:));
    

else
    
    if ndims(u)==3
        ns = size( u , 3);
        ub  = squeeze(u(:,:,end));
        u  = squeeze(sum( u , 3 ) / ns);                
    else
        
        ub  = squeeze(u);
    end
    ub = ub';
    ub = ub(msk,:);
    u = u';
    U   = u;
    eta = eta';
    D   = eta +  ones(nt,1) * d;
    u   = u(msk,:);
    U   = U(msk,:);
    eta = eta(msk,:);
    D   = D(msk,:);
    depmsk = D < options.drylim;
    

end

% Flooding and drying catch

ustar = squeeze( ustar);
ustar = ustar';
ustar = ustar(msk,:);


ustar(depmsk) = 0;
tau = -sign(ub).*ustar.^2;




%Calcu
d_dx.eta = derivative( eta , diff(x) );
d_dx.U   = derivative( U   , diff(x) );
d_dx.u   = derivative( u   , diff(x) );

d_dx.D   = derivative( D   , diff(x) );

dd = D -eta;
Qd = u.*D ./ dd;
d_dx.Qd   = derivative( Qd   , diff(x) );

if options.includeVerticalGradientTerms
   %
   sigma = -((1 : 1 : ns) - 0.5)/ns;
   sigma =repmat( sigma , [size(eta,1) , 1 ,size(eta,3) ] );
   d_dx.s = sigma.* d_dx.D + d_dx.eta;

   %
else
    d_dz = [];
end

end

%

% function dFuncds = vertDerivative( func , var)
% % Function to calculate vertical gradients
%   %  nt = size(func,1);
%   %  ns = size(func,2);
%   %  nx = size(func,3);
%     
%   if strcmpi( var , 'p' )
%         % Variable is p-like.
%         
%   elseif strcmpi( var , 'u' )
%         % Variable is u-like.
%         
%   elseif strcmpi( var , 'w' )
%         % Variable is w-like.
%         
%   end  
% 
% end

%
%
%
function dFuncdx = derivative( func , dx )
    %
    nt = size(func,1);

    if ndims(func) == 3
        ns = size(func,2);
        nx = size(func,3);
    else
        ns = 1;
        nx = size(func,2);
    end
            
    if ns == 1
        dFuncdx = zeros( nt,nx );
        %
        for jj=2:size(func,2)-1
            %  
            dFuncdx(:,jj) = ( func( :,jj+1 ) - func(:,jj-1) ) ./ (dx(jj-1) +dx(jj));
            %
        end
        %
        dFuncdx(:,1) = nan;
        dFuncdx(:,end) = nan;
    else
        dFuncdx = zeros( nt,ns,nx );
        for jj=2:size(func,3)-1
            %        
            dFuncdx(:,:,jj) = ( func( :,:,jj+1 ) - func(:,:,jj-1) ) ./ (dx(jj-1) +dx(jj));
            %
        end
        %
        dFuncdx(:,:,1) = nan;
        dFuncdx(:,:,end) = nan;
    end
end

function [f , C] = cospec(  x , y  , df, fs  )
%
nd = ndims(x);
%
if nd == 3
    nt     = size(x,1);
    nx     = size(x,3);
    nsigma = size(x,2);
    x      = reshape( x , [nt , nx * nsigma ] );
    y      = reshape( y , [nt , nx * nsigma ] );
end
%
options.overlap = 0.5;
options.nfft = min( floor( fs / df ) , size( x , 1 ) );


window = hamming( options.nfft );    
overlap = floor( options.overlap * options.nfft );



%remove nans
x(isnan(x)) = 0;
y(isnan(y)) = 0;

%calculate cospectra
[ C , f ] = cpsd( x , y , window , overlap , options.nfft , ...
    fs , 'onesided');

if nd == 3
    %Depth integrate if vertically distributed
    nf = length(f);
    C  = reshape( C , [nf , nsigma, nx ] );
    C  = squeeze(sum(C , 2) / nsigma);
end
C = real(C);
end

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