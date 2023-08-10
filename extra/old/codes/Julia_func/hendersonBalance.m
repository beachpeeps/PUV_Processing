function [f , F , dF_dx , S , E] = hendersonBalance( t , x , u ,  eta , d , tstart , tend , varargin )
%
options.type = 'NH';
options.df   = 0.01;
options.fs   = 1/mean(diff(t));
options.g    = 9.81;
options.DepthAveraged = true;


options.type = 'NLSWE';

options = parseOptions( options , varargin );

d(d<0.1) = nan;

%
[t,u,U,eta, d,D,d_dx , Qd ] = procVars( t,x,u,eta,d,tstart,tend,options);

[f , F ,  S , E] = NLSWEBalance( x ,U , eta  , d , D ,Qd, d_dx ,options);




%Flux Gradient
names = fieldnames(F);
for jj=1:length(names)
    name = names{jj};
    dF_dx.(name) = derivative( F.(name) , diff(x) );
end

%
end

function [f , F ,  S , E]  = NLSWEBalance(  x, u, eta  , d , D ,Qd, d_dx, options)
%ShortHands
g  = options.g;
df = options.df;
fs = options.fs;

%
% cospectra - used in energy
%
if ndims(eta)==3
    eta0 = squeeze(eta(:,1,:));
else
    eta0 = eta;
end

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

%keyboard


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
S.sw   = dQddxSxx;
S.tot  = S.sw;
%

end

function [t,u,U,eta,d,D, d_dx , Qd ] = procVars( t,x,u,eta,d,tstart,tend,options)
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
    u   = permute(u , [2,3,1]);
    eta = eta';
    D   = eta + ones(nt,1) * d;
    
    eta = repmat(eta , [1,1,ns] );
    eta = permute(eta , [1,3,2] );
    
    D = repmat(D , [1,1,ns] );
    D = permute(D , [1,3,2] );

else
    
    if ndims(u)==3
        ns = size( u , 3);
        u  = squeeze(sum( u , 3 ) / ns);
    end
    u = u';
    U   = u;
    eta = eta';
    D   = eta +  ones(nt,1) * d;
    U   = U(msk,:);
    eta = eta(msk,:);
    D   = D(msk,:);
end

dd = D -eta;
Qd = U.*D ./ dd;
Qd(dd<0.01) = nan;





%Calcu
d_dx.eta = derivative( eta , diff(x) );
d_dx.U   = derivative( U   , diff(x) );
d_dx.u   = derivative( u   , diff(x) );
d_dx.Qd   = derivative( Qd   , diff(x) );

d_dx.D   = derivative( D   , diff(x) );
end



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