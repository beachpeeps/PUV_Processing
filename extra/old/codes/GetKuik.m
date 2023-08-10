
function [mdir1,mdir2,spr1,spr2,skw,kur]=GetKuik(a1,b1,a2,b2)


% calculation of "estimator free" directional parameters from 
%  Kuik, A. J., G. Ph. van Vledder, and L.H. Holthuijsen, A method
%  for routine analysis of pitch-and-roll buoy data. J. Phys.
%  Oceanogr., 18, 1020-1034, 1988.

% first moment mean wave direction
mdir1=atan2(b1,a1)*(180/pi);

% a1 b1 spread
spr1=2*(1-sqrt(a1.^2+b1.^2));
spr1=sqrt(spr1)*180/pi;

% --- 2nd moment mean direction and spread -----

% a1b1 mean direction, in radians
mdir1=atan2(b1,a1);
% kuik et. al. a2-b2 second moment m2, relative to a1b1 mean dir
m2=a2.*cos(2*mdir1)+b2.*sin(2*mdir1);

% a2b2 mean direction
mdir2=0.5*atan2(b2,a2)*(180/pi);

% convert a1b1 mean direction from radians to degrees
mdir1=mdir1*(180/pi);

% a2b2 mean dir has 180 deg amiguity. find one that is closest to
% a1b1 mean dir.
tdif=abs(mdir1-mdir2);
mdir2(tdif > 90)=mdir2(tdif > 90)-180;
mdir2(mdir2 < 0)=mdir2(mdir2 < 0)+360;

% calculate a2-b2 spread after kuik
spr2=sqrt((1.0-m2)/2)*(180/pi);

%----  turn mdir 1 negative directions in positive directions
%       before returning
mdir1(mdir1 < 0)=mdir1(mdir1 < 0)+360;

%--- skewness and kurtosis calcs 

% m1 after Kuik et. al.
rm1=sqrt(a1.^2+b1.^2);
% 2 times mean direction
t2=2*atan2(b1,a1);
% n2 after Kuik et. al.
rn2=b2.*cos(t2)-a2.*sin(t2);
% m2 after Kuik et. al.
rm2=a2.*cos(t2)+b2.*sin(t2);

% kurtosis
kur=(6.-8.*rm1+2.*rm2)./((2*(1.-rm1)).^2);
% skewness
skw=-rn2./(.5*(1-rm2)).^1.5;

end
