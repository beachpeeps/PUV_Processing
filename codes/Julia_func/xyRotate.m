function varargout = xyRotate(x,y,theta,xo,yo)
% modified from Nathaniel Plant
% clockwise rotation around origin XO,YO
% xyRotate  Rotate data.
%   [XR YR] = xyRotate(X,Y,THETA) rotates the coordinates X,Y by THETA
%   (theta in radians). (Origin does not move, X0 = 0, Y0 = 0)
%
%   [XR YR] = xyRotate(X,Y,THETA,XO,YO) rotates the coordinates around the
%   origin XO,YO.

%%
% Default origin is 0,0.
if exist('xo','var')==0; xo=0; end
if exist('yo','var')==0; yo=0; end

% Rotation Matrix
A = [cos(theta) sin(theta); -sin(theta) cos(theta)];

% Translate and rotate.
out = [x(:)-xo y(:)-yo]*A;

if nargout==1
   varargout = {out};
elseif nargout==2
   varargout(1) = {reshape(out(:,1),size(x))};
   varargout(2) = {reshape(out(:,2),size(y))};
end
