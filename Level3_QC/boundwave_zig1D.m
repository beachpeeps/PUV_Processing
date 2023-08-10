%% Boundwave timeseries generation - can be changed to 2D
% require x (pressure timeseries) input
%%
function [Z_IG] = boundwave_zig1D(x, depth)
    %% Implimenting SWASH boundwave BC for boundwave - 1D
    
    N = length(x);
    fs = 2;
    freq_full = [0:fs/N:fs/2 -fs/2+fs/N:fs/N:-fs/N];

    
    khmin = 0.265; % corresponds to f = 0.04Hz
    khmax = 10; % corresponds to f=0.48Hz
    eta = zeros(N,1);
    phase3_store = zeros(N,N);
    A3_store = zeros(N,N);
    l_store = NaN(N,N);
    
    d = depth;
    g = 9.81;
    A = fft(x,N);
    amp = abs(A)/N; 
    amp(2:end-1) = 2*amp(2:end-1);
    phase = angle(A); 
    df = fs/N;
    lcheck = [];
    
    for idf1 = 1:N % compute boundwave
        A1 = amp(idf1);
        f1 = freq_full(idf1);
        phase1 = phase(idf1);
        omega1 = 2*pi*f1;
        k1 = get_wavenumber(omega1,d);
        
        if k1*d > khmin && k1*d < khmax && A1 ~=0
    %  *********************   first loop over f1             XXXXXXXXXXXXXXXXXXXXXXXXXXXZ
            for idf2 = idf1+1:N
                A2 = amp(idf2);
                f2 = freq_full(idf2);
                phase2 = phase(idf2);
                omega2 = 2*pi*f2;
                k2 = get_wavenumber(omega2,d); 
                if k2*d > khmin && k2*d < khmax && A2 ~=0
    %%%%%%%%%%%%%%   2nd loop over f2 %%%%%%%%%%%%%%%%%% f2 and f3 are changing.  
                    f3 = abs(f2 - f1);
                    cosdt1 = -1;
                    k3 = real(sqrt(k1.^2+k2.^2+2.*k1.*k2.*cosdt1));
                    phase3 = phase2 - phase1 + pi;                  % 100% bound  ++++++++++++++++++
                    omega2 = -omega2;  
                    cg3 = 2*pi*f3/k3;
                    C = (omega1+omega2) * ( (omega1*omega2).^2 /(g^2) - k1*k2*cosdt1 )...
                            - 0.5 * ( omega1*k2.^2/(cosh(k2*d).^2) + omega2*k1.^2/(cosh(k1*d).^2) );
                    T1 = -g*k1*k2*cosdt1 ./ (2.*omega1*omega2); % equation 11, first half
                    T2 = 1/2/g.*(omega1.^2+omega2^2+omega1.*omega2); %equation 11, 2nd half
                       
                    T3 = g*(omega1+omega2) ./...
                            (( g*k3*tanh(k3*d) - (omega1+omega2).^2 ).*(omega1*omega2));
                        
                    D1D = T1+T2+T3*C; 
    
                    % SWASH code stuff 
                    A3 = abs(D1D).*A1.*A2;                     %%%%    *****  A3 is real and positive  Predicted Bound wave amp for a single SS pair *****
                    l = round(f3/df)+1;                      %%%%  l is ig band number - need to allow for 0 frequency in FFT array
                    eta(l) = eta(l) + A3.*exp(sqrt(-1).*phase3);       %%%%%  eta(l) is  COMPLEX ....sum-over all f2, f1 with f3= f2 -f1.
                    phase3_store(idf1, idf2) = phase3;
                    A3_store(idf1, idf2) = A3;
                    l_store(idf1, idf2) = l;
                end
            end % idf2
        end
       % loop over f1 where f3 = f2 - f1
    end % idf1
    eta_backup = eta;
    
    %% Get timeseries
    eta_norm = zeros(1, N);
    fSS = find(abs(freq_full) > 0.04 | abs(freq_full) < 0.004);
    eta_norm = eta*(N/2);
    eta_norm(fSS)=0;
    Z_IG = ifft(eta_norm,'symmetric');
end
