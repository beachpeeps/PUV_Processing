%% Boundwave timeseries generation 
% require x (eta incoming timeseries) input

%%
function [Z_IG, E_bound, f_bound, alpha] = boundwave_zig2D(x, depth, FC, fm, data_input)
    %% Implimenting SWASH boundwave BC for boundwave - 2D
    N = length(x);
    fs = 2;
    freq_full = [0:fs/N:fs/2 -fs/2+fs/N:fs/N:-fs/N];
    
    khmin = 0.265; % corresponds to f = 0.04Hz
    khmax = 10; % corresponds to f=0.48Hz
    eta = zeros(N,1); phase3_store = zeros(N,N); A3_store = zeros(N,N); l_store = NaN(N,N);
    
    d = depth;
    g = 9.81;
    A = fft(x,N);
    aa = abs(A).^2 / N;
    amp(1) = aa(1);
    amp(2:N/2+1) = 2*aa(2:N/2+1);
    phase = angle(A); 
    df = fs/N;
    lcheck = [];

    %% get direcional spectra
    dd(:,:)= mem_est(FC.a1, FC.a2, FC.b1, FC.b2);
    for ii = 1:size(dd,2)
        dd_interp(:,ii)=interp1(fm, dd(:,ii), freq_full);
    end
    for i=1:length(amp) % loop through freq bands
        ds(i,:)=dd_interp(i,:)*amp(i); % mutiply by the freq band total energy
    end  
    if strcmp(char(data_input),'PUV')
        ds = [ds(:,90:end) ds(:,1:89)];
    elseif strcmp(char(data_input),'MOP')
        ds = fliplr(ds);
    end
            
    clear e1 e2 e12
    omega = 2*pi*freq_full;
    kwav = get_wavenumber(omega,d);
    
    %cosine of difference angle
    theta1 = [180:360]; 
    theta2= theta1; lt = length(theta1);
    [THETA1,THETA2] = meshgrid(theta1,theta2);
    cosdt = cosd(THETA2-THETA1+180);
    
    
    iidf = floor(0.004/df):1:ceil(0.04/df);% changed from 10 to 1 interval
    f_bound = df*(iidf);
    i_swell = find(freq_full > 0.04 & freq_full < 0.25);         
%     figure
%     polarPcolor(freq_full(1:i_swell(end)),[1:360],ds(1:i_swell(end),:)')
%%
    clear E_bound eta
    for idf = iidf % index deltaf, 1st loop runs through everything on the df
       % idf
            E3 = NaN(i_swell(end)-idf, size(cosdt,1), size(cosdt,1)); E3_denorm = E3;
            E0 = E3; E1 = E3; E2 = E3;
            phase3 = NaN(i_swell(end)-idf,1); A3 = NaN(i_swell(end)-idf,1);
            parfor if1 = i_swell(1):i_swell(end)-idf % for all swell indices, 2nd loop runs through everything on f1, so all combos of df,f1 and f2 = f1+df
                % f1
                f1 = freq_full(if1);
                phase1 = phase(if1);
                omega1 = omega(if1);
                k1 = kwav(if1);

                % f2
                f2 = freq_full(if1 + idf);
                phase2 = phase(if1 + idf);
                omega2 = -omega(if1 + idf);
                k2 = kwav(if1 + idf);

                % get bound wave wavenumber k
                k3 = sqrt(k1.^2+k2.^2+2*k1.*k2.*cosdt);
                C = (omega1+omega2) * ( (omega1*omega2).^2 /(g^2) - k1*k2*cosdt ) - 0.5 * ( omega1*k2.^2/(cosh(k2*d).^2) + omega2*k1.^2/(cosh(k1*d).^2) );
                T1 = -g*k1*k2*cosdt ./ (2.*omega1*omega2); % equation 11, first half
                T2 = 1/2/g.*(omega1.^2+omega2^2+omega1.*omega2); %equation 11, 2nd half
                T3 = g*(omega1+omega2) ./ (( g*k3.*tanh(k3.*d) - (omega1+omega2).^2 ).*(omega1*omega2));
                D = T1+T2+T3.*C;
            
                % get energy at f2,theta1 and f1,theta2
                e1 = ds(if1,theta1);
                e2 = ds(if1+idf,theta2);
                % get energy at f3, boundwave, need squared values because energy
                E0(if1,:,:) = (e1'*e2);
                E1(if1,:,:) = D.^1.*(e1'*e2);
                E2(if1,:,:) = D.^2.*(e1'*e2);

                E3(if1,:,:) = D.^2.*(e1'*e2); % *****  A3 is real and positive  Predicted Bound wave amp for a single SS pair *****
                
                phase3(if1) = phase2 - phase1 + pi;
                if phase3(if1) > pi; phase3(if1) = phase3(if1) - 2*pi; end

            end % f1 loop

            E_bound(idf) = 2*nansum(E3*df, 'all');
            num1 = nansum(E1, [2 3]); num1 = nansum(num1(idf:end));
            num2 = nansum(E2, [2 3]); num2 = nansum(num2(idf:end));
            denom = nansum(E0, [2 3]); denom = nansum(denom(idf:end));
            M1(idf) = num1 ./ denom;
            M2(idf) = num2 ./ denom;
            A3 = nansum(E3,[2 3]); % integrate over dtheta1 dtheta2
            A3 = A3 * (N/2)^2;
            A3 = sqrt(A3);

            eta(idf) = nansum(A3.*exp(sqrt(-1).*phase3).*df); % integrate over df      %%%%%  eta(df) is  COMPLEX ....sum-over all f2, f1 with f3= f2 -f1.
          
    end
    alpha = M2 ./ (M1.^2);

    f_bound = freq_full(1:iidf(end));
    %% Get timeseries
    fSS = find(abs(freq_full) > 0.04 | abs(freq_full) < 0.004);
    eta_norm = zeros(1, length(freq_full));
    eta_norm(iidf) = eta(iidf);
    eta_norm(fSS)=0;
    
    Z_IG = ifft(eta_norm,'symmetric');
    
end
