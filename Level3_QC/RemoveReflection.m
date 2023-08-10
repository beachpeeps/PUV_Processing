function [Zp,Zr,Ztrend,Znl] = RemoveReflection( Z, U, H,  dt, varargin )
    %
    % Remove reflection component (prop. in negative x direction) from a data signal to obtain only the
    % forward propagating component (in the positive x direction). Method is
    % essentially a non-linear and dispersive correctio to the original
    % Herbers & Guza methodology. 
    %
    % Z = Z(+) + Z(-)   ( Z progressive + Z reflected)
    % U = U(+) + U(-)   ( U "           + U "        )
    % 
    
    %There are three options:
    % 
    % 1) Herbers & Guza
    %-------------------------------------
    %
    % This method is based on linear wave theory and estimates the 
    % reflected and incident systems as
    %
    % Zp = .5 Z + .5 sqrt(H/g) U (Eq. 1)
    % Zr = .5 Z - .5 sqrt(H/g) U (Eq. 2)
    %
    % 2) Sheremet extension
    %-------------------------------------
    %
    % In this case we add full linear dispersion so that in Fourier space
    % we have
    % 
    % Z+ = .5 * Z + .5 * U / F
    % Z- = .5 * Z - .5 * U / F
    %
    % The factor F accounts for the depth attenuation of the signal and is 
    % given as
    %
    % F = cosh( kH ) / cosh( kh + kz )
    %
    % Alternatively, if U is the depth averaged velocity it is given as
    %
    % F = k D / w
    %
    % 3) Nonlinear weakly dispersive theory
    %--------------------------------------
    % We only consider this for low-frequency components (i.e. f<fmax), so we
    % consider (Eq. 1) in the frequency domain (denoted by _h)
    %
    % Zp_h = .5 Z_h + .5 sqrt(H/g) U_h (from Eq. 1)   For f<fmax
    %                                                                   (Eq. 3)                                                                3)
    % Zp_h = Z_h                                      else
    % 
    % and Z = ifft(Zp);
    %
    %--------------------------------------------------------------------------
    
    
    % 
    % INPUT PARSING/Init
    %
	p      = inputParser;
    
    grav = 9.81;
    addParameter( p , 'z' , -H );
    addParameter( p , 'method' , 'HG' );
    addParameter( p , 'depav' , false );
    addParameter( p , 'g' , grav );
    addParameter( p , 'flim' , [-1,-1] );
    addParameter( p , 'trend' , true );
    p.parse( varargin{:} );
    
    z      = p.Results.z;
    method = p.Results.method;
    depav  = p.Results.depav;    
    g      = p.Results.g;
    flim   = p.Results.flim;
    trend   = p.Results.trend;
               

    %% Detrend signals
    
    if trend
        %
        Ztrend = Z-detrend(Z);
        Z      = Z - Ztrend;
        U      = detrend( U );
        %
    else
        %
        Ztrend = 0 * Z;
        %
    end

    % Frequency range
    % NOTE: We only correct the positive frequencies (Below Fnyq = 1/(2dt)). For
    % the negative frequencies we assume that Z was a real signal, such that
    % they are complex conjugates.

    
    nt     = size(Z,1);
    nx     = size(Z,2);
    T      = nt*dt;
    df     = 1/T;
    f      = (0:1:nt-1)*df;    
    fnyq  = floor( (nt+1)/2 )*df;
    ifnyq = floor( (nt+1)/2 );
    
    % Bandfiltering
    fmax = flim(2);
    fmin = flim(1);
    
    %
    if (fmax >= fnyq || fmax < 0 )
        %
        imax = ifnyq;
        %
    else
        %
        [dum,imax] = min(abs(f-fmax));
        %
    end
    
    if (fmin  < 0 )
        %
        imin = 2;
        %
    else
        %
        [dum,imin] = min(abs(f-fmin));
        %
    end    
    
    fmsk  = imin:imax;

    %Calculate Transforms
    Zh   = fft(Z);
    Uh   = fft(U);       
    Zph = Zh;
    Zrh = 0*Zh;

    if length(H)==1
        %
        D = H*ones(nt,nx);
        %
    else
        %
        if size(H,1) == 1
            %
            D = H*ones(nt,1);
            %
        elseif size(H,2) == 1
            %
            D = H'*ones(nt,1);
            %
        else
            %
            D = H;
            %
        end            
        %
    end
    
    w    = 2*pi*f;
    fact = zeros( size( Zh ) );
    nl = false;
    switch lower(method)
        %
        case {'hg'}
            %
            % Herbers and guza
            %
            fact = sqrt(D/g);
            %
        case {'sher'}
            %
            % Sheremet method
            %
            if depav 
                %
                % Depth averaged velocity as input
                %
                for ix = 1 : nx
                    %
                    k = disper(w,D(1,ix));

                    fact(:,ix) = (k.*D(1,ix))./w;
                    fact(1,ix) =0;                
                    %
                end
                %
            else
                %
                % Local velocity
                %                
                for ix = 1 : nx
                    %
                    k = disper(w,D(1,ix));
                    d = D(1,ix);
                    fact(:,ix) = w ./ g ./ k .* cosh( k.*d ) ./ cosh( k.*d + k.*z );
                    fact(1,ix) =0;                
                    %
                end  
                %
            end
            %
        case {'smit'}
            %
            nl = true;
            Zhnl = fft(Z.^2);
            factnl = zeros( size( Zh ) );
            %
            if depav
                %
                for ix = 1 : nx
                    %                
                    d = D(1,ix);

                    alpha = - d ./ g .* ( 1. / 6. );

                    fact(:,ix) = sqrt( d./g ) ./ ( 1 + w.^2 .* alpha );
                    fact(1,ix) =0;
                    factnl(:,ix) =1 ./ ( 8.*d .* ( 1 + w.^2 .* alpha ) );
                    factnl(1,ix) =0;                
                    %
                end
                %
            else
                %
                for ix = 1 : nx
                    %                
                    d = D(1,ix);

                    alpha = d ./ g .* ( 1. / 6. + z./d.^2 .* ( z./2 + d ) );

                    fact(:,ix) = sqrt( d./g ) ./ ( 1 + w.^2 .* alpha );
                    fact(1,ix) =0;
                    factnl(:,ix) =1 ./ ( 8.*d .* ( 1 + w.^2 .* alpha ) );
                    factnl(1,ix) =0;                
                    %
                end
                %
            end
            %
    end
    %   
    Zph(fmsk,:) = Zh(fmsk,:)/2 + fact(fmsk,:).*Uh(fmsk,:)/2;
    Zrh(fmsk,:) = Zh(fmsk,:)/2 - fact(fmsk,:).*Uh(fmsk,:)/2;
    %
    Znl = real(0*Zph);
    %
    if nl
        %
        
        Znl(fmsk,:) = factnl(fmsk,:).*Zhnl(fmsk,:);
        Znl = ifft(Znl,'symmetric');
        
        Zph(fmsk,:) = Zph(fmsk,:) + factnl(fmsk,:).*Zhnl(fmsk,:);
        Zrh(fmsk,:) = Zrh(fmsk,:) - factnl(fmsk,:).*Zhnl(fmsk,:);
        %
    end
    %
    % Inverse transform
    %
    Zp = ifft(Zph,'symmetric');
    Zr = ifft(Zrh,'symmetric');
    %
end




function k=disper(w, h, g)
    % DISPER   Linear dispersion relation
    %
    % usage:   k  = disper(w,h,g)     or
    %
    %          k  = wave number             (2 * pi / wave length)
    %          w  = wave angular frequency  (2 * pi / wave period)
    %          h  = water depth
    %          g  = gravitational acceleration constant, optional (DEFAULT 9.81)
    %
    %          relative error in k*h < 2.5e-16 for all k*h
    %

    %          programmer: G. Klopman, Delft Hydraulics, 6 Dec 1994

    if nargin < 3,
      g = 9.81;
    end;

    w2 = max(0.0001,(w.^2) .* h ./ g);
    q  = w2 ./ (1 - exp (-(w2.^(5/4)))) .^ (2/5);

    for j=1:2,
      thq     = tanh(q);
      thq2    = 1 - thq.^2;
      a       = (1 - q .* thq) .* thq2;
      b       = thq + q .* thq2;
      c       = q .* thq - w2;
      arg     = zeros (size (q));
      iq      = find (a ~= 0);
      arg(iq) = (b(iq).^2) - 4 .* a(iq) .* c(iq);
      arg(iq) = (-b(iq) + sqrt(arg(iq))) ./ (2 * a(iq));
      iq      = find (abs(a.*c) < 1.0e-8 * (b.^2));
      arg(iq) = - c(iq) ./ b(iq);
      q       = q + arg;
    end;

    k = sign(w) .* q ./ h;

    ik    = isnan (k);
    k(ik) = zeros(size(k(ik)));
end

