
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




%Calcu
d_dx.eta = derivative( eta , diff(x) );
d_dx.U   = derivative( U   , diff(x) );
d_dx.u   = derivative( u   , diff(x) );
d_dx.Qd   = derivative( Qd   , diff(x) );

d_dx.D   = derivative( D   , diff(x) );
end