function [m2,m3] = sectormoments( T, sig , flim ,inorm)
    %----------------------------------------------------------------------
    % Description
    %----------------------------------------------------------------------
    % A simple script to calculate the second and third-order "sector"
    % moments from a given wave signal. For the second-order statistics we
    % define the "sector-moment" merely as the variance contained within a
    % frequency band. For the third order moment, the sector-moments are
    % basically integrals over sub-regions of the bi-spectrum (regions as 
    % in de Bakker et al. 2015, JPO).
    %
    %----------------------------------------------------------------------
    % VERSION Author   Date        Comment
    %----------------------------------------------------------------------
    % 1       P.B.Smit 25-Feb-2017 Initial release
    % 2       J Fiedler 21-Sep-2018 Modifications to norms
    %
    %----------------------------------------------------------------------
    % Input: 
    %----------------------------------------------------------------------
    % T  : [1,1] REAL
    %      lenght/duration of the signal
    % sig: [nt,nx] REAL
    %      "surface elevation signal" (or other scalar property)
    %      consiting of nt sampled observations in time at a regular
    %      interval. For nx>1 the sectormoments will be calculated 
    %      for each column of z seperately.
    % flim: [1 : NSEC+1 ] REAL
    %       vector consisting of nlim frequencies [f1,f2,..., fnlim] which
    %       define the sector boundaries ( e.g. frequencies between 
    %       f1 < f < f2 are considered the first sector, f2 < f < f3 the
    %       second sector etc.
    % inorm [1,1]: INTEGER OPTIONAL
    %       (optional): chooses the normalization for the third order
    %        moments. 
    %        [=0]: (DEFAULT) use the second-order moment of the total signal
    %        [=1]: Normalize with relevant second-order moments (see note)
    %----------------------------------------------------------------------
    % ORDER OF CALCULATION
    %---------------------------------------------------------------------- 
    % 1) Detrend Signals
    % 2) Perform the passband filters
    % 3) Perform Hilbert transform on the filter signal
    % 4) Calculate m2
    % 5) Calculate m3
    % 6) Check if the results "add up" (skewness from input signal is sum
    %    of skewness in m3
    % 7) Return output
    %
    %----------------------------------------------------------------------
    %                           IMPLEMENTATION
    %----------------------------------------------------------------------
    if nargin < 4
        %
        inorm = 0;
        %
    end
    
    if nargin >3
        %disp('hello!!')
        inorm;
    end
    %
    % 1 Remove trends...
    %----------------------------------------------------------------------
    sig = detrend(sig);

    %
    % 2 Hilbert Transform
    %----------------------------------------------------------------------
    sigf = hilbert( sig );    
    
    %
    % 3 Passband filter 
    %----------------------------------------------------------------------
    sigf = filter( T,sigf,flim);
        
    %
    % 4 Calculate the second order sector moments
    %----------------------------------------------------------------------
    m2 = mom2( sigf );
    %
    % 5 Calculate the third order sector moments
    %----------------------------------------------------------------------
   [ m3 ] = mom3( sigf , inorm , m2, flim,T);
    
   
   
   
    %
    % 6 A simple check on the calculations if the total skewness and
    % variance from the sector calculations adds up to the
    % skewness/variance directly form the signal
    %----------------------------------------------------------------------
    Var = mean( sig.^2,1)';
    Sk  = mean( sig.^3,1)' ./ Var.^(3/2);
    
    DifVar =  max(abs(1-abs(Var./m2(:,  end))));
    DifSk  =  (abs(1-abs(real(m3.tot')./Sk)));
    %
    % Note for "tiny" signals, normalized skewness very unreliable anyway - ignore
    %
    DifSk( abs(Sk) < 0.001 ) = 0;
    DifSk = max(DifSk);    
    %
    if ( flim(1) <= 0. && (flim(end)/ (1/T) > floor(size( sig,1)/2) + 1) )
        %
        %
        % Only check if all the sectors combined comprise the whole signal
        % else there is no point...
        if (DifVar > 0.01)
            %
            disp('Warning, sums do not add up in variance...' );
            %
        end
        %
        if (DifSk > 0.05)
            %
            disp('Warning, sums do not add up in skewness...' );
            %
        end    
        %
    end
    %
    % 7 Return
    %----------------------------------------------------------------------
    % Done
end
%==========================================================================
%--------------------------------------------------------------------------
%                               HELPER FUNCTIONS
%--------------------------------------------------------------------------
%==========================================================================
function [ m3 ] = mom3( zf , inorm , var , flim,TT)
    %
    nx     = size(zf,2);
    nl     = size(zf,3)-1;
    %
    %
    npairs = (nl - 1)*nl / 2;

    
    m3.self  = zeros( nx , nl );
    m3.sum   = zeros( nx , npairs );
    m3.dif   = zeros( nx , npairs );
    m3.norm3   = zeros( nx , npairs );
    

    pairs = zeros( 2 , npairs );
    
    ipair = 0;
    for ii = 1 : nl - 1
        %        
        for jj = ii +1  : nl
            %
            ipair = ipair + 1;
            pairs( : , ipair ) = [ ii , jj ];
            %
        end
        %
    end
    
        
    %
    % Total number of interactions triplets - for each triplet 1!
    % 
       
    itrip = 0;
    for ii = 1 : nl - 2
        %        
        for jj = ii +1  : nl - 1
            %
            for kk = jj +1  : nl
                %            
                itrip = itrip + 1;               
                %
            end
            %
        end
        %
    end    
    
       
    ntrip = itrip;
    m3.trip  = zeros( nx , ntrip );
    
    itrip = 0;
    trip = zeros( 3,ntrip);
    for ii = 1 : nl - 2
        %        
        for jj = ii +1  : nl - 1
            %
            for kk = jj + 1  : nl
                %            
                itrip = itrip + 1;
                trip( : , itrip ) = [ ii , jj , kk];
                %
            end
            %
        end
        %
    end     
    
    
    %%
    % Factors for the doubles
    %
    fac = [3/2; ...
           3/4; ...
           3/4; ...
           3/4];
       
    if ( inorm == 3)
        %
        fac(1:4) = 1;
        %
    end       
       
    %
    % For each sum/difference pair do:
    %    
    tot = 0;
    for ipair = 1 : npairs
        %
        
        ip1 = pairs( 1 , ipair );
        ip2 = pairs( 2 , ipair );
        
        comb = [ ip1,ip2,ip2 ; ...  % ig - ss - ss
                 ip1,ip1,ip2];      % ss - ss - ss                    
        %
        for jj = 1 : 2
            %
            T = mean( zf( : , : , comb(jj,1) ) .* zf( : , : , comb(jj,2) ) .* conj( zf( : , : , comb(jj,3) ) ) );
            tot = tot + T*fac(jj);                        
            
            switch inorm
                %
                % Normalizations...
                %
                case {1}
                    %
                    norm = (sqrt( var( :,comb(jj,1) ).*var( :,comb(jj,2) ).*var( :,comb(jj,3) ) ))';
%                     if jj==1
%                     disp(T)
%                     disp(norm)
%                     end
                    %
                case {2}
                    %
                    norm = 1;
                    %
                case {3,4}
                    %
                    % The new normalizations
                    %
                    if jj == 2
                        %
                        % For sum interactions this is weird
                        %
                        sig1=      zf( : , : , comb(jj,1) );                    
                        sig2=      zf( : , : , comb(jj,1) ).*conj( zf( : , : , comb(jj,3) ));
                        
                        %f= filter( TT,(sig2),flim(ip1:ip1+1) );
                        z_tmpr = filter( TT,real(sig2),flim(ip1:ip1+1) );
                        z_tmpr = z_tmpr + conj(z_tmpr);
                        z_tmpi = filter( TT,imag(sig2),flim(ip1:ip1+1) );
                        z_tmpi = z_tmpi + conj(z_tmpi);                       
                        z_tmp = z_tmpr + 1i * z_tmpi;
                        z_tmp = z_tmp-mean(z_tmp);
                        %keyboard
                        if inorm == 3
                            %
                            sig2 = z_tmp(:,:,1);
                            %
                        end
                        
                        n2 = mean( sig2.*conj(sig2) );
                        n1 = mean( sig1.*conj(sig1) );
                        %
                    else
                        %
                        % Difference interactions
                        %
                        sig1=      zf( : , : , comb(jj,1) );                    
                        sig2=      zf( : , : , comb(jj,3) ).*conj( zf( : , : , comb(jj,3) ));
                        sig2 = sig2-mean(sig2);
                        z_tmp = filter( TT,sig2,flim(ip1:ip1+1) );

                        z_tmp = z_tmp + conj(z_tmp);
                        %
                        if inorm == 3
                            %                            
                            sig2 = z_tmp(:,:,1); 
                            sig2 = sig2-mean(sig2);
                            %
                        end

                        %
                    end
                    n1 = mean( sig2.*conj(sig2) );
                    n2 = 0.5*mean( sig1.*conj(sig1) );
                    %                                       
                    norm = sqrt(n1 .* n2);
                    m3.norm3(:,ipair) = norm;
                otherwise
                    %
                    % Default norm.
                    %
                    norm = (sqrt( var( :,end ).*var( :,end ).*var( :,end ) )');
            end

            if jj == 2
                m3.sum(:,ipair) = fac(jj) * T ./ norm;
            else
                %disp([fac(jj) ipair])
              
                m3.dif(:,ipair) = fac(jj) * T ./ norm;
            end
            
            %
        end
        %
    end
    %
    % The "self" interactions
    %
    for ipair = 1 : nl
        %
        
        ip1 = ipair;
        ip2 = ipair;
        
        comb = [ ip1,ip2,ip2 ; ...  % ig - ss - ss
                 ip1,ip1,ip2];      % ss - ss - ss                    
        %
        %
        T = 3/4*mean( zf( : , : , ip1 ) .* zf( : , : , ip1 ) .* conj( zf( : , : , ip1 ) ) );

        tot = tot + T;                        

        switch inorm
            case {1}
                norm = (sqrt( var( :,comb(jj,1) ).*var( :,comb(jj,2) ).*var( :,comb(jj,3) ) ))';
            case {2}
                norm = 1; 
%             case {3}
%                 %                
%                 sig1 =zf( : , : , ip1 )*3/4;
%                 sig2 = conj(zf( : , : , ip1 )) .* ( zf( : , : , ip1 ) );               
%                 sig2 = sig2 - mean(sig2);
%                 
%                 z_tmp = filter( TT,sig2,flim(ip1:ip1+1) );                
%                 z_tmp = z_tmp + conj(z_tmp);                
%                 
%                 sig2 = z_tmp(:,1);
% 
%                 n1 = mean(sig1.*conj(sig1));
%                 n2 = mean(sig2.*conj(sig2));               
%                 norm = sqrt(n1*n2)/sqrt(2);
            
            otherwise
                norm = (sqrt( var( :,end ).*var( :,end ).*var( :,end ) )');
        end
        %
        m3.self(:,ipair) = T ./ norm;
        %
    end    
    
    
    %
    % For each triplet do
    %
    for itrip = 1 : ntrip
        %        
        ip1 = trip( 1 , itrip );
        ip2 = trip( 2 , itrip );
        ip3 = trip( 3 , itrip );
        %        
        T = mean( zf( : , : , ip1 ) .* zf( : , : , ip2 ) .* conj( zf( : , : , ip3 ) ) );
        tot = tot + T*3/2;
        %
        switch inorm
            %
            case {1}
                %
                norm = (sqrt( var( :,ip1 ).*var( :,ip2 ).*var( :,ip3 ) ))';
                %
            case {2}
                norm = 1;                
            otherwise
                %
                norm = (sqrt( var( :,end ).*var( :,end ).*var( :,end ) )');
                %
        end
        %
        m3.trip(:,itrip) = 3/2 * T ./ norm;
        %
    end    
    %
    tot = tot ./ (var( : , end).^(3/2))';    
    
    m3.tot = tot;
    m3.itrip = trip;
    m3.ipair = pairs;
    %       
end

function m2 = mom2( zf )
    %
    nt     = size(zf,1);
    nx     = size(zf,2);
    nl     = size(zf,3);    
    %        
    m2 = zeros( nx , nl );
    for jj = 1 : nl
       %
       % Calculate the variances for normalizing
       %
       m2( : , jj ) = .5 * mean( zf( : , : , jj) .* conj( zf( : , : ,jj ) ) ,1 )';
       %
    end    
    %
end

function zf = filter( T,z,flim)
    %
    nx = size( z, 2 );
    nt = size( z, 1 );
    df = 1 / ( T );
    nl = length(flim) - 1;    
    %
    zf = zeros( nt , nx , nl + 1 );
    %
    z = z - ones(nt,1)*mean(z,1);
    FourierAmp = fft( z );
    %
    for ilim = 1 : nl
       %
%        jf = floor( flim(ilim : ilim+1)/df )+1;
        dt = T/nt;
        fnyq = 1/2/dt;
        lim = flim(ilim:ilim+1);
        lim(2) = min(fnyq,lim(2));
        jf = floor(lim/df)+1;
       %       
       tmp = FourierAmp;
       tmp( 1:jf(1) , :)   = 0;
       tmp( jf(2):end,: ) = 0;
       
       zf( : , : , ilim ) = ifft( tmp );
       %
    end
    zf( : , : , nl+1) = sum( zf( : , : , 1:nl) , 3 );
    %
end


