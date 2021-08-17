clear
load('../data/temp.mat','PUV_process');

%%



for i=1:length(PUV_process)
try
Hsig(i) = PUV_process(i).Hsig.Hs;
fmean(i) = PUV_process(i).ids.fcentroid_swell;
end
end



%%



id = 883;
fm = PUV_process(id).Spec.fm;
i_swell = PUV_process(id).ids.i_swell;
%%
FC = PUV_process(id).FC;
df = fm(i_swell(2))-fm(i_swell(1));
g = 9.81;
d = PUV_process(id).Eflux.depth;
%% get MEM estimate

clear dd ds
%compute directional spreding function using MEM estimator

    
    dd(:,:)= mem_est(FC.a1, FC.a2, FC.b1, FC.b2);

    for i=1:length(fm) % loop through freq bands
        ds(i,:)=dd(i,:)*PUV_process(id).Spec.SSE(i); % mutiply by the freq band total energy
    end

%% do the bound wave thang
dds = repmat(ds,1,1,360); % double the angles, double the fun
clear e1 e2 e12
omega = 2*pi*fm;
kwav = get_wavenumber(omega,d);


theta1 = [1:30 330:360]; theta2 = theta1;
 
e1D = PUV_process(id).Spec.SSE(1:length(fm));
[THETA1,THETA2] = meshgrid(theta1,theta2);

dtheta = pi/180;

iidf = 1:20:round(0.1/df);
f_bound=df*(iidf);
E_bound = zeros(1,max(iidf));

%cosine of difference angle
cosdt = cosd(THETA2-THETA1+180);

 tic   
 hw = waitbar(0./length(f_bound),'Looping through all IG difference frequencies');
 
 for idf = iidf %index deltaf, 1st loop runs through everything on the df
     waitbar(idf./max(iidf),hw,['df = ' num2str(fm(idf),'%2.4f') ' Hz'])
%         n = 0;
     for ii = i_swell(1):i_swell(end)-idf % for all swell indices, 2nd loop runs through everything on f1, so all combos of df,f1 and f2 = f1+df
%          n = n+1;
         omega2 = omega(ii+idf); % m
         k2 = kwav(ii+idf);
         
         omega1 = -omega(ii); % n
         k1 = kwav(ii); %note that I do not define k1 based on omega here, 
                        %because it is only a sign switch on omega!!
         
         % get bound wave wavenumber k
         k3 = sqrt(k1.^2+k2.^2+2*k1.*k2.*cosdt);
       
        %interaction coefficient party time
           %Following Hasselmann 
        C = 1i*(omega1+omega2) * ( (omega1*omega2).^2 /(g^2) - k1*k2*cosdt )...
            - 1i*0.5 * ( omega1*k2.^2/(cosh(k2*d).^2) + omega2*k1.^2/(cosh(k1*d).^2) );
        T1 = -g*k1*k2*cosdt ./ (2.*omega1*omega2); % equation 11, first half
        T2 = 1/2/g.*(omega1.^2+omega2^2+omega1.*omega2); %equation 11, 2nd half
        
        T3 = g*(omega1+omega2) ./...
            (( g*k3*tanh(k3*d) - (omega1+omega2).^2 ).*(omega1*omega2));
        
        
        D = abs(T1+T2+T3*C);    
        
%         1D version (I am too lazy to write good code).
        cosdt = -1;
            % get bound wave wavenumber k
         k3 = sqrt(k1.^2+k2.^2+2*k1.*k2.*cosdt);
        C = 1i*(omega1+omega2) * ( (omega1*omega2).^2 /(g^2) - k1*k2*cosdt )...
            - 1i*0.5 * ( omega1*k2.^2/(cosh(k2*d).^2) + omega2*k1.^2/(cosh(k1*d).^2) );
        T1 = -g*k1*k2*cosdt ./ (2.*omega1*omega2); % equation 11, first half
        T2 = 1/2/g.*(omega1.^2+omega2^2+omega1.*omega2); %equation 11, 2nd half
        
        T3 = g*(omega1+omega2) ./...
            (( g*k3*tanh(k3*d) - (omega1+omega2).^2 ).*(omega1*omega2));
        
        
        D1D = abs(T1+T2+T3*C); 
        
     % get energy at f2,theta1 and f1,theta2
     e1 = squeeze(dds(ii+idf,theta1,theta2)); % E(f2,all-thetas)
     e2 = squeeze(dds(ii,theta1,theta2)); % E(f1,all-thetas) 
     
     % get energy at f3, boundwave, need squared values because energy
     E3(ii,:,:) = D.^2.*e2.*e1;
     E31D(ii) = D1D.^2.*e1D(ii).*e1D(ii+idf);
       
     end %swell loop
       %using only 1 dtheta here, because the first summation over all angles
       %is normalized. I think.
       E_bound(idf) = 2*sum(E3(2:end,:,:)*df*dtheta,'all');
       E_bound1D(idf) = 2*sum(E31D*df,'all');

     clear E3 k3
     
 end % df loop
 toc
 close(hw)
E_bound = E_bound(iidf);
E_bound1D = E_bound1D(iidf);
%
%%
clf
esector = sum(dds(:,theta2,theta1),[2 3]).*dtheta;

semilogy(fm,esector)
hold on

semilogy(fm,sum(ds,2))
hold on
semilogy(f_bound,E_bound)
semilogy(f_bound,E_bound1D)

 xlabel('f (Hz)')
 ylabel('E (m^2/Hz)')
 title(['id = ' num2str(id) ' and not a date because I am lazy']) 
 legend('sum of E(f,theta_{calc})','sum of E(f,theta) TOTAL','2D bound','1D bound')
 saveas(gcf,['../viz/test_boundwave_' num2str(id) '.png'])

toc