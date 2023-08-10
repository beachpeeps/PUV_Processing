
function zf = filter_zeta( T,z,flim)
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
       jf = floor( flim(ilim : ilim+1)/df )+1;
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


