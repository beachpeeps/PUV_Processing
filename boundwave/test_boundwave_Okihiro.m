clear
% load('temp.mat','PUV_process');
load('/Users/juliafiedler/Downloads/IB_data.mat')

fm = double(Fq_IB);
df = fm(2)-fm(1);
iswell = find(fm<=0.25);
% i_swell = find(fm>0.04 & fm<0.25);
%
%
% for id = 1:length(a1_IB);
%     Hs(id) = 4*sum(Ed_IB(:,id).*df);
%     a1 = a1_IB(:,id);
%     a2 = a2_IB(:,id);
%     b1 = b1_IB(:,id);
%     b2 = b2_IB(:,id);
%     try
%         mdir1(:,id)=atan2(b1,a1)*(180/pi); %mean dir 1st mom
%         % 1st mom spread
%         spr1_temp=2*(1-sqrt(a1.^2+b1.^2));
%         spr1(:,id)=sqrt(spr1_temp)*180/pi;
%     end
% end
% %
id = 4385;
% id = 300;
spec1D = Ed_IB(:,id);
% a1 = a1_IB(:,id);
% a2 = a2_IB(:,id);
% b1 = b1_IB(:,id);
% b2 = b2_IB(:,id);
% % %
% %
% kk = isnan(a1);
% a1(kk) = [];
% b1(kk) = [];
% a2(kk) = [];
% b2(kk) = [];
% fmk = fm(kk);
% spec1D
%add noise
% n = length(fm);
% a1 = a1+ 0.2*randn(n,1);
% a2 = a2+ 0.2*randn(n,1);
% b1 = b1+ 0.2*randn(n,1);
% b2 = b2+ 0.2*randn(n,1);



% %
% clear dd ds
% %compute directional spreading function using MEM estimator
%
%     dd(:,:)= mem_est(a1, a2, b1, b2);
%     for i=1:length(fmk) % loop through freq bands
%         ds(i,:)=dd(i,:)*spec1D(i); % mutiply by the freq band total energy
%     end
%
% pcolor(fm,1:360,ds'); shading flat
%% do the bound wave thang
% dds = repmat(ds,1,1,360); % double the angles, double the fun

% id = 323;
% fm = PUV_process(id).Spec.fm;
% i_swell = PUV_process(id).ids.i_swell;

%%
df = Bw_IB(1);
g = 9.81;
d = depth_IB;


Hs = 4*sum(Ed_IB.*df);
%% do the bound wave thang
omega = 2*pi*fm;
kwav = get_wavenumber(omega,d);

clear idf

iidf = 1:round(0.15/df);
f_bound=df*(iidf);
E_bound = zeros(1,max(iidf));

% dtheta  = 1/360; 
dtheta = 1;
%cosine of difference angle

% just putting theta in here for shiggles
theta1 = 42; 
theta2 = 42;
lt1 = length(theta1);

[THETA1,THETA2] = meshgrid(theta1,theta2);
cosdt = cosd(THETA2-THETA1+180);

tic
hw = waitbar(0./length(f_bound),'Looping through all IG difference frequencies');

for idf = iidf %index deltaf, 1st loop runs through everything on the df
    waitbar(idf./max(iidf),hw,['df = ' num2str(fm(idf),'%2.4f') ' Hz'])
    for ii = 1:length(fm)-idf % for all swell indices, 2nd loop runs through everything on f1, so all combos of df,f1 and f2 = f1+df
        
        % Here we are computing D(f+f',-f',\deltatheta+180), following
        % Herbers 94 where f is the difference-frequency
        omega2 = omega(ii+idf); % this is f+f'
        k2 = kwav(ii+idf);

        
        omega1 = -omega(ii); % this is -f'
        k1 = kwav(ii);
        


        
        k3 = sqrt(k1.^2+k2.^2 + 2*k1*k2*cosdt);
%         k3 = sqrt(k1.^2+k2.^2 +k1*k2*cosdt); %Herbers k (is wrong?)

        %Following Hasselmann 
        C = 1i*(omega1+omega2) * ( (omega1*omega2).^2 /(g^2) - k1*k2*cosdt )...
            - 1i*0.5 * ( omega1*k2.^2/(cosh(k2*d).^2) + omega2*k1.^2/(cosh(k1*d).^2) );
        T1 = -g*k1*k2*cosdt ./ (2.*omega1*omega2); % equation 11, first half
        T2 = 1/2/g.*(omega1.^2+omega2^2+omega1.*omega2); %equation 11, 2nd half
        
        T3 = g*(omega1+omega2) ./...
            (( g*k3*tanh(k3*d) - (omega1+omega2).^2 ).*(omega1*omega2));
        
        
        D = abs(T1+T2+T3*C);  
        
        % get energy at f2 and f1
        %          e1 = squeeze(dds(ii+idf,theta1,theta2)); % E(f2,all-thetas)
        %         e2 = squeeze(dds(ii,theta1,theta2)); % E(f1,all-thetas)
        e1 = spec1D(ii);
        e2 = spec1D(ii+idf);
        
        % get energy at f3, boundwave, need squared values because energy
        %     	E3(ii,:,:) = D.^2*e2.*e1;
        E3(ii) = D.^2*e1.*e2;
        

    end %swell loop

    E_bound(idf) = 2*sum(E3*df);%*df*dtheta);

end % df loop

toc
close(hw)
E_bound = E_bound(iidf);
%
% for ii = i_swell(1):i_swell(end)-1
%     pcolor(squeeze(THETA3(ii,:,:))); shading flat
%     title(num2str(ii))
%     caxis([-3 3])
%     pause(0.1)
% end
%

% plot check
%
%
clf
% esector = sum(dds(:,theta2,theta1)*dtheta,[2 3]);

semilogy(fm,spec1D)
hold on
% semilogy(fm,esector)


semilogy(f_bound,E_bound)
xlim([0 0.25])

xlabel('f (Hz)')
ylabel('E (m^2/Hz)')
title(['id = ' num2str(id) ' and not a date because I am lazy'])
legend('spec1D','calculated 1D bound wave energy')
saveas(gcf,['test_1Dboundwave_' num2str(id) '.png'])
%

%%
%  Plot C
% xUnits = omega(1:end-1).*sqrt(d/g);
%
% loglog(xUnits,abs(D*d))