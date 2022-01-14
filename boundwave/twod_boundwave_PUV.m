clear
cd '/Users/athinalange/Desktop/IG_BC_remote_work/PUV_Processing-julia'
%load('temp.mat','PUV_process');
addpath('1stQC')
addpath('codes')
addpath('boundwave')
%%



for i=1:length(PUV_process)
try
Hsig(i) = PUV_process(i).Hsig.Hs;
fmean(i) = PUV_process(i).ids.fcentroid_swell;
spread(i) = PUV_process(i).dir.spread2_ss_sum;
end
end

%%

%id = 795%1582; %1581;
ids = [123 256 1143 1126 1580 274 795 376];%1580; %1581;
%figure(8);clf
for jj = 1:8%1:length(PUV_process)
%clearvars -except boundwave PUV_process Hsig fmean spread jj eta_level ids
id = ids(jj)
if isempty(PUV_process(id).Spec) == 0
sprintf('hi')
fm = PUV_process(id).Spec.fm;
i_swell = PUV_process(id).ids.i_swell;
%%
FC = PUV_process(id).FC;
df = fm(i_swell(2))-fm(i_swell(1));
g = 9.81;
d = PUV_process(id).Eflux.depth;
%% get MEM estimate

clear dd ds T1 T2 T3 D1D D E3 E31D E_bound E_bound1D e1D dds ds dd
%compute directional spreding function using MEM estimator

    dd(:,:)= mem_est(FC.a1, FC.a2, FC.b1, FC.b2);

    for i=1:length(fm) % loop through freq bands
        ds(i,:)=dd(i,:)*PUV_process(id).Spec.Spp(i); % mutiply by the freq band total energy
    end    

% shift dd to work in 'coming from' angle because I am backwards today
ds = [ds(:,90:end) ds(:,1:89)];

%  Fake '2D' version   
% ds = zeros(901, 360);
% ds(:, 270) = PUV_process(id).Spec.SSE;
% figure(4);clf
% pcolor(fm, [1:360],ds')
% shading flat
%% do the bound wave thang
dds = repmat(ds,1,1,360); % double the angles, double the fun
clear e1 e2 e12
omega = 2*pi*fm;
kwav = get_wavenumber(omega,d);

theta1 = [180:360]; 
%theta1 = [270-29:270+30];
%theta1 = [1:360];% theta2 = theta1; lt = length(theta1);
theta2= theta1; lt = length(theta1);

e1D = PUV_process(id).Spec.Spp(1:length(fm));
[THETA1,THETA2] = meshgrid(theta1,theta2);

dtheta = 1; %not used

iidf = 1:10:round(0.1/df);
f_bound=df*(iidf);
% E_bound = zeros(1,max(iidf));

%cosine of difference angle
cosdt = cosd(THETA2-THETA1+180);

 tic   
 hw = waitbar(0./length(f_bound),'Looping through all IG difference frequencies');
 
 for idf = iidf %index deltaf, 1st loop runs through everything on the df
     idf
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
        C = (omega1+omega2) * ( (omega1*omega2).^2 /(g^2) - k1*k2*cosdt )...
            - 0.5 * ( omega1*k2.^2/(cosh(k2*d).^2) + omega2*k1.^2/(cosh(k1*d).^2) );
        T1 = -g*k1*k2*cosdt ./ (2.*omega1*omega2); % equation 11, first half
        T2 = 1/2/g.*(omega1.^2+omega2^2+omega1.*omega2); %equation 11, 2nd half
        
        T3 = g*(omega1+omega2) ./...
            (( g*k3.*tanh(k3.*d) - (omega1+omega2).^2 ).*(omega1*omega2));
        
        
        D = T1+T2+T3.*C; % aa(idf,ii) = D(270,270);   
        
%         1D version (I am too lazy to write good code).
        cosdt1 = -1;
            % get bound wave wavenumber k
         k3 = sqrt(k1.^2+k2.^2+2*k1.*k2.*cosdt1);
        C = (omega1+omega2) * ( (omega1*omega2).^2 /(g^2) - k1*k2*cosdt1 )...
            - 0.5 * ( omega1*k2.^2/(cosh(k2*d).^2) + omega2*k1.^2/(cosh(k1*d).^2) );
        T1 = -g*k1*k2*cosdt1 ./ (2.*omega1*omega2); % equation 11, first half
        T2 = 1/2/g.*(omega1.^2+omega2^2+omega1.*omega2); %equation 11, 2nd half
        
        T3 = g*(omega1+omega2) ./...
            (( g*k3*tanh(k3*d) - (omega1+omega2).^2 ).*(omega1*omega2));
        
        
        D1D = T1+T2+T3*C; %aa1d(idf,ii) = D1D; 
        
     % get energy at f2,theta1 and f1,theta2
%      e1 = squeeze(dds(ii+idf,theta1,theta2)); % E(f2,all-thetas)
     e1 = ds(ii+idf,theta2);
%      e2 = squeeze(dds(ii,theta1,theta2)); % E(f1,all-thetas) 
     e2 = ds(ii,theta1);
     
     % get energy at f3, boundwave, need squared values because energy
     E3(ii,:,:) = D.^2.*(e1'*e2);
     E31D(ii) = D1D.^2.*e1D(ii).*e1D(ii+idf);
       
     end %swell loop
       %adding a 2pi here because i like it better this way
       E_bound(idf) = 2*sum(E3*df,'all');
       E_bound1D(idf) = 2*sum(E31D*df,'all');

    
 end % df loop
 toc
 close(hw)
E_bound = E_bound(iidf);
E_bound1D = E_bound1D(iidf);

boundwave(id).E_bound = E_bound;
boundwave(id).E_bound1D = E_bound1D;
boundwave(id).f_bound = f_bound;
boundwave(id).ds = ds;
boundwave(id).fm = fm;
%
%%
% 
% esector = 1/lt.*sum(dds(:,theta2,theta1),[2 3]);
% subplot(6,2,1+(jj-1)*2)
% semilogy(fm,esector,'linewidth',4)
% hold on
% semilogy(fm,sum(ds,2),'linewidth',2)
% hold on
% semilogy(f_bound,E_bound,'linewidth',4)
% semilogy(f_bound,E_bound1D,'linewidth',2)
% 
%  xlabel('f (Hz)')
%  ylabel('E (m^2/Hz)')
%  title(['id = ' num2str(id) ', Hs = ' num2str(Hsig(id),'%2.2f') ' m, spread = ' num2str(spread(id),'%2.2f') '\circ']) 
%  legend('sum of E(f,theta_{calc})','sum of E(f,theta) TOTAL','2D bound','1D bound')
%  set(gca, 'FontSize', 15)
% ax2 = subplot(6,2,2+(jj-1)*2)
% % polarPcolor(R,theta,Z,'Ncircles',3)
% 
% [h,c] = polarPcolor(fm',0:360,[ds ds(:,end) ],'Nspokes',9,'typeRose','meteo'); shading flat
% hold on
% c.Label.String  = 'E(f,\theta)'; c.Label.FontSize = 12;
% % polarplot(theta1(1)*ones(length(fm)),fm)
% % theta1 = [1:15 345:360];
% dr = ds;
% dr(i_swell,theta1) = 1;
% set(gca, 'FontSize', 15)
% [h2,c2] = polarPcolor(fm',0:360, [dr dr(:,end)],'Nspokes',9,'typeRose','meteo','colbar',0); shading flat
% h2.FaceAlpha = 0.2;
% %  xlabel('f (Hz)')
% %  ylabel('Dir (\circ), 0\circ  = onshore')
% h =  title('E(f,\theta), 270\circ shorenormal'); h.Position(2) = 1.1;

% 
% figure(2);clf
%  esector = 1/lt.*sum(dds(:,theta2,theta1),[2 3]);
% subplot(1,2,1)
% semilogy(fm,esector,'linewidth',4)
% hold on
% semilogy(fm,sum(ds,2),'linewidth',2)
% hold on
% semilogy(f_bound,E_bound,'linewidth',4)
% semilogy(f_bound,E_bound1D,'linewidth',2)
% xlabel('f (Hz)')
% ylabel('E (m^2/Hz)')
% legend('sum of E(f,theta_{calc})','sum of E(f,theta) TOTAL','2D bound','1D bound')
% title(['id = ' num2str(id) ', Hs = ' num2str(Hsig(id),'%2.2f') ' m, spread = ' num2str(spread(id),'%2.2f') '\circ']) 
% 
% set(gca, 'FontSize', 15)
% ax2 = subplot(1,2,2)
% % polarPcolor(R,theta,Z,'Ncircles',3)
% [h,c] = polarPcolor(fm',0:360,[ds ds(:,end) ],'Nspokes',9,'typeRose','meteo'); shading flat
% hold on
% c.Label.String  = 'E(f,\theta)'; c.Label.FontSize = 12;
% % polarplot(theta1(1)*ones(length(fm)),fm)
% % theta1 = [1:15 345:360];
% dr = ds;
% dr(i_swell,theta1) = 1;
% set(gca, 'FontSize', 15)
% [h2,c2] = polarPcolor(fm',0:360, [dr dr(:,end)],'Nspokes',9,'typeRose','meteo','colbar',0); shading flat
% h2.FaceAlpha = 0.2;
% %  xlabel('f (Hz)')
% %  ylabel('Dir (\circ), 0\circ  = onshore')
% h =  title('E(f,\theta), 270\circ shorenormal'); h.Position(2) = 1.1;
% % % saveas(gcf,['../viz/test_boundwave_' num2str(id) '.png'])
end
    
toc
end
