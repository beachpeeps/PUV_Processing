function Ur = get_ursellNumber(f,h,H)
% This function gets the Ursell number given:
% f: frequency
% h: water depth
% H: wave height
% Output: Ur = a/h/(kh)^2

omega = 2*pi*f;
k = get_wavenumber(omega,h);

a = H/2;
Ur = (a./h)./(k.*h).^2;
