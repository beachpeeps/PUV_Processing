clear
load('forBill')
a1_og=a1;
a2_og=a2;
b1_og=b1;
b2_og=b2;

%% cos2s
clear a1 a2 b1 b2
a1 = a1_og+0.2*rand(length(a1_og),1);
a2 = a2_og+0.2*rand(length(a2_og),1);
b1 = b1_og+0.2*rand(length(b1_og),1);
b2 = b2_og+0.2*rand(length(b2_og),1);

kk=isnan(a1);
fm(kk)=[];
a1(kk)=[];
b1(kk)=[];
a2(kk)=[];
b2(kk)=[];
spec1D(kk)=[];

%compute directional spreding function using MEM estimator
dd(:,:)= mem_est(a1, a2, b1,b2);

    for i=1:length(fm) % loop through freq bands
        ds(i,:)=dd(i,:)*spec1D(i); % mutiply by the freq band total energy
    end

% shift dd to work in 'coming from' angle because I am backwards today
ds = [ds(:,90:end) ds(:,1:89)];
figure(3);clf
pcolor(fm, [1:360], ds')
shading flat
%%
mdir1=atan2d(b1,a1);
%sigma = rad2deg(sqrt(2*(1 - sqrt(a1.^2 + b1.^2))));
%s = 2./(sigma.^2) - 1;
%
% calculate s based on eq 9.22 at
% https://www.sciencedirect.com/topics/engineering/directional-spreading
r1=sqrt(a1.^2 + b1.^2);
s=(r1*pi)./(1-r1*pi);
% cos2sC
for ii = 1:length(fm)
    % keep angles between -180 and  +180 (+/-90 when divided by 2)
    ang=[1:360] - mdir1(ii);
    ang(ang > 180)=ang(ang > 180)-360;ang(ang < -180)=ang(ang < -180)+360;
    C = 2.^(2*s(ii)).*(gamma(s(ii)+1)).^2./(2*pi*gamma(2*s(ii)+1));
    cos2s(ii,:) = C*(cosd(ang/2).^(2*s(ii))); 
    ds(ii,:) = cos2s(ii,:) * spec1D(ii);
end
ds = [ds(:,90:end) ds(:,1:89)];

figure(1);clf
pcolor(fm, [1:360], cos2s') % Directional Distribution
shading flat

figure(2);clf
pcolor(fm, [1:360], ds')
shading flat