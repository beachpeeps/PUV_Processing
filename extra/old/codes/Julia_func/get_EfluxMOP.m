function [EfluxX, EfluxY] = get_EfluxMOP(MOP)
% [EfluxX, EfluxY] = get_EfluxMOP(MOP)
% this function calculates energy flux from the MOP line, given the input
% structure MOP, containing frequency array and its bounds, depth, 1D
% spectra, fourier coeffs a1 and b1, and the shorenormal angle.


% define constants
rho = 1028;    % density (kg/m^3)
g = 9.81; % gravity (m/s^2)

%% Estimate Energy Flux
% spectral estimator from T.H.C. Herbers, using group velocity
k = get_k(MOP.frequency,MOP.depth);
Cg = get_cg(k,MOP.depth);   %group speed, 0.5 * omega.* k .* ( 1 + (2*k*depth)./sinh(2*k*depth) );


% these use the fourier coefficients from a North=0 deg grid
posX= rho * g * Cg' .* MOP.a1 .* MOP.spec1D;
posY = rho * g * Cg' .* MOP.b1 .* MOP.spec1D;
% units = kg/m^3 * m/s^2 * m/s * m^2 *s = kg*m/s^2

% convert to radians and rotate into +x = 0 on shorenormal
theta = deg2rad(-MOP.shorenormal);
% the following does clockwise rotation
[posXR, posYR] = xyRotate(posX,posY,theta);

% sum over all frequencies (is this right??)
df = double(MOP.fbounds(2,:)-MOP.fbounds(1,:))'; % f bounds are not constant!!
EfluxX = sum(posXR' .* df);
EfluxY = sum(posYR' .* df);


%% helper functions

    function varargout = xyRotate(x,y,theta,xo,yo)
        % modified from Nathaniel Plant
        % clockwise rotation around origin XO,YO
        % xyRotate  Rotate data.
        %   [XR YR] = xyRotate(X,Y,THETA) rotates the coordinates X,Y by THETA
        %   (theta in radians). (Origin does not move, X0 = 0, Y0 = 0)
        
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
    end

    function k = get_k(f,h)
        
        %k = getk(f,h)
        %Credit F. Fedderson
        % returns the wavenumber of the gravity wave
        % dispersion relation, by using newtons method
        % the initial guess will be the shallow water wavenumber
        
        omega = 2*pi*f;
        k = omega./sqrt(g*h);
        f = g*k.*tanh(k.*h) - omega.^2;
        
        while max(abs(f))>1e-10
            
            dfdk = g*k.*h.*(sech(k.*h)).^2 + g*tanh(k.*h);
            k = k - f./dfdk;
            f = g*k.*tanh(k.*h) - omega.^2;
        end
    end



    function cg = get_cg(k,h)
        % Falk Feddersen (c) 2001
        %
        % function that takes the vector wavenumbers and depths
        % and returns the group velocity
        %
        % function cg = funwaveC_get_cg(k,h)
        % check to make sure k & h are the same size
        [nk,mk] = size(k);
        [nh,mh] = size(h);
        
        if (nh==1)
            if (mh==1)
                cg = 0.5*(g*tanh(k*h)+g*(k.*h).*(sech(k.*h).^2))./sqrt(g*k.*tanh(k.*h));
                return;
            end
        end
        
        if (nk ~= nh)
            disp('** Error k & h have wrong sizes');
            return;
        end
        
        if (mk ~= mh)
            disp('** Error k & h have wrong sizes');
            return;
        end
        
        cg = 0.5*(g*tanh(k.*h)+g*k.*h.*(sech(k.*h).^2))./sqrt(g*k.*tanh(k.*h));
    end
end
