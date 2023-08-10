function MOPname = get_MOPname(lat)
% MOPname = get_MOPname(lat,lon)
% retrieves the sitelabels (e.g. 'D0666', 'D0045') for given latitude input
% by using a nearest neighbor search. This might produce unexpected results
% if the input is wildly off from the actual MOP lat. Probably should
% write a catch statement for this.
% Note also that the latitude needs to have 4 decimal places, as the MOP
% lines are separated by approx
% 
% INPUT: 
% lat, lon : in decimals, can be a vector
%
% OUTPUT:
% MOPname: string array specifying the site label
%
% Example:
% MOPname = get_MOPname(32.5735, -117.1491)
% will result in MOPname = 'D0045';
%
% Julia Fiedler, 2019 jfiedler@ucsd.edu

% get all the lat/lon/site info for the MOPS lines in SoCal
url = 'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/model/MOP_alongshore/socal_alongshore_nowcast.nc';
Latitude = ncread(url,'metaLatitude'); 
% Longitude = ncread(url,'metaLongitude'); 
SiteLabel = ncread(url,'metaSiteLabel');

indLat = knnsearch(Latitude(:),lat(:));

% TODO: Get Longitude to work as well??
% indLon = knnsearch(Longitude(:),lon(:));
% indMOP = intersect(indLat,indLon);


for i=1:length(indLat)
    MOPname(i) = deblank(convertCharsToStrings(SiteLabel(:,indLat(i))));
end

% MOPname = deblank(convertCharsToStrings(SiteLabel(:,indMOP)));




