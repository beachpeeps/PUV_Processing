function [xCor,normCor] = crossCorr2d(Nlag,mat1,x2)
% Nlag = 1001; %make odd for going over 0 zero lag.
% Nloc = 850;

[Nt,Nloc] = size(mat1);
xCor = nan(Nlag,Nloc);

for ix = 1:Nloc %for all x locations
    for ilag = 1:Nlag % we want lags to go -500:500 
        ii = Nlag-ilag+1:Nt-ilag+1; % ii goes from 1001-(1:1001):6003-(1:1000)
        jj = ilag:Nt-Nlag+ilag; % jj goes from (1:1001):6002-1001+(1:1001)
        
        xCor (ilag, ix) = mean(mat1(ii,ix) .*x2(jj));
    end
end

% norm for 
% indices change for normalization bc lags make it shorter
L = length(x2);
inds = 1+floor(Nlag/2):L-floor(Nlag/2);
n1 = mean(mat1(inds,:).*conj(mat1(inds,:)));
n2 = mean(x2(inds).*conj(x2(inds)));
normCor = n1.*n2;
