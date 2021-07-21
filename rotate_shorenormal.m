%% Rotate to shorenormal coordinates
%
% From ENU coordinate system, pulling shorenormal angle from CDIP,
% rotating to onshore, shorenormal (+x) and alongshore, north (+y) [Right-handed coordinate system]
%
% Input: MOP: double
%        U,V: velocities in buoy coordinates [1 x time] 
%
%% Author: 
% Athina Lange, SIO July 2021


function [Uprime, Vprime, shorenormal] = rotate_shorenormal(U, V);

load('moplist') % fulll list of all CDIP buoy's as of 07/22/19
[indx,tf] = listdlg('ListString',moplist, 'SelectionMode','single', 'InitialValue',[1], 'Name', 'What MOP do you want to use?');
stn = char(moplist(indx));
ncfile = strcat('http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/model/MOP_alongshore/', string(stn(end-4:end)), '_nowcast.nc');
shorenormal=ncread(ncfile,'metaShoreNormal'); % mop shore normal

% Coordinate system is in ENU, so rotate +x to North, then CCW rotation to
% MOP shorenormal, then rotate 180 to get 'going to' frame of reference. 
rotation = rem(90+(360 - shorenormal)+180, 360); % +x is Onshore, Shorenormal, +y is North, alongshore
alpha_rad = rotation*pi/180; 
rot_mat   = [cos(alpha_rad)  -sin(alpha_rad);
             sin(alpha_rad)  cos(alpha_rad)];

uv=rot_mat*[U';V'];
Uprime=uv(1,:).'; 
Vprime=uv(2,:).';


end
