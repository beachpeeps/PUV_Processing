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
options.overlap = 0;
options.nfft = 7200;%min( floor( fs / df ) , size( x , 1 ) );


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