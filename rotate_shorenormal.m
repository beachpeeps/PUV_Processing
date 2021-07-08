%% Rotate to shorenormal coordinates
%
% From a buoy coordinate system (WNU), pulling shorenormal angle from CDIP,
% rotating to onshore, shorenormal (+x) and alongshore, south (+y)
%
% Input: MOP: double
%        U,V: velocities in buoy coordinates [1 x time] 
%


function [Uprime, Vprime, shorenormal] = rotate_shorenormal(MOP, U, V);

url = ['https://cdip.ucsd.edu/mops/?moplist=San_Diego_County&pub=public&mop=D0' char(string(MOP)) '&xitem=info']
name = 'Normal';
shorenormal = urlfilter(url,name)

rotation = shorenormal - 90; % +x is Onshore, Shorenormal, +y is South, alongshore
alpha_rad = rotation*pi/180; 
rot_mat   = [cos(alpha_rad)  sin(alpha_rad);
            -sin(alpha_rad)  cos(alpha_rad)];

uv=rot_mat*[U';V'];
Uprime=uv(1,:).'; 
Vprime=uv(2,:).';


end