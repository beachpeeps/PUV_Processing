function [A,nens] = calculate_fft2(X,nfft,fs)
        
        [n,m]=size(X);
        
        num = floor(2*n/nfft)-1;
        nens = num;
        %[n nfft num]
        
        X = X-ones(n,1)*mean(X); %demean
        
        X(isnan(X)) = 0;
        X = detrend(X);
        X = X-ones(n,1)*mean(X); %demean
        
        sumXt = (X'*X)/n; %get fourier coeffs
        
        %WIN = hanningwindow(@hamming,nfft);
        jj = [0:nfft-1]';
        WIN = 0.5 * ( 1 - cos(2*pi*jj/(nfft-1)));
        
        
        A = zeros(num,nfft);
        
        % set it up so that SQR(|A(i,:)|^2) = sum(X^2) (Parseval's Thm)
        
        varXwindtot = 0;
        
        for i=1:num,
            istart = (i-1)*(nfft/2)+1;
            istop = istart+nfft-1;
            Xwind = X(istart:istop);
            Xwind = Xwind - mean(Xwind);  % demean.   detrend?
            Xwind = detrend(Xwind); %detrend
            lenX = length(Xwind); %why is this not evident from nfft
            varXwind =( Xwind'*Xwind)/lenX; %get variance, normalized
            varXwindtot = varXwindtot + varXwind; %add to total variance
            Xwind = Xwind .* WIN; %window it
            tmp = ( Xwind'*Xwind)/lenX; % get windowed covariance
            if (tmp == 0),
                Xwind = Xwind * 0.0;
            else
                Xwind = Xwind*sqrt(varXwind/tmp); %parseval's thm
            end;
            A(i,:) = fft(Xwind')/sqrt(nfft);
            meanA2 = mean( A(i,:) .* conj(A(i,:)));
        end;
    end