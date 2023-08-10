function Ho = reverseshoal(H,h,f)

k = get_wavenumber(2*pi*f,h);
kh = k .*h;

% shoaling coefficient, 
% Nielsen, P. (2009) Coastal and Estuarine Processes. Vol. 29, Advanced
% Series on Ocean Engineering, World Scientific, 343 pp.
Ks = 1./ sqrt( tanh(kh).*(1.+2*kh./sinh(2*kh)) );


Ho = H./Ks;