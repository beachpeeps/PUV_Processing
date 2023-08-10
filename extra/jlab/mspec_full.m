function[varargout]=mspec_full(varargin)

if strcmp(varargin{1},'--t')
    mspec_test;return
end
if strcmp(varargin{1},'--f')
    mspec_fig;return
end

if isscalar(varargin{1})
  deltat=varargin{1};
  varargin=varargin(2:end);
else
  deltat=1;
end
  

%Sort out input arguments
lambda=1;  %This means use the average multitaper spectrum
if isstr(varargin{end})
    if strfind(varargin{end},'ada')
        lambda=varargin{end-1};
        varargin=varargin(1:end-2);
    end
end

x=varargin{1};
psi=varargin{end};
na=length(varargin);

if na==2
    if isreal(x)
        y=[];
    else
        y=conj(x);
    end
elseif na==3
    y=varargin{2};
end

varargout=mspec_one(deltat,x,y,psi,lambda);


    
function[cellout]=mspec_one(dt,x,y,psi,lambda)

%In real cases, multiply by sqrt(2) to get one-sided spectrum
if isempty(y)        %One real-valued
     [f,mmatx]=mtrans(sqrt(2)*x,psi);
else
     if ~isreal(x)   %Two complex-valued
        [f,mmatx,mmaty]=mtrans(x,y,psi);
     else            %Two real-valued
        [f,mmatx,mmaty]=mtrans(sqrt(2)*x,sqrt(2)*y,psi); 
     end
end


if isempty(y) %One time series
     if lambda==1
         cellout{2}=avgspec(mmatx,mmatx).*dt;
     else
         var=squared(vstd(x,1)); 
         cellout{2}=adaptspec(abs(mmatx).^2,lambda,var).*dt;
     end
else         %Two time series
     if lambda==1
        cellout{2}=avgspec(mmatx,mmatx).*dt;
        cellout{3}=avgspec(mmaty,mmaty).*dt;
        cellout{4}=avgspec(mmatx,mmaty).*dt;
     else
        %For two time series one should do the adaptive spectra on both
        %with the same coefficients
        var=squared(vstd(x,1))+squared(vstd(y,1)); 
        
        [s,dk]=adaptspec(abs(mmatx).^2+abs(mmaty).^2,lambda,var);
        cellout{2}=squeeze(frac(sum(dk.^2.*abs(mmatx).^2.*dt,2),sum(abs(dk).^2,2)));
        cellout{3}=squeeze(frac(sum(dk.^2.*abs(mmaty).^2.*dt,2),sum(abs(dk).^2,2)));
        cellout{4}=squeeze(frac(sum(dk.^2.*mmatx.*conj(mmaty).*dt,2),sum(abs(dk).^2,2))); 
     end
end


%Corrections for zero component, and for Nyquist with even and odd length
%Both zero are Nyquist are shared for even length time series, so divide both by two.
%Only zero is shared for odd length, since the Nyquist does not appear. 
%This is best visualized by drawing N equally spaced points on the unit circle.

for i=2:length(cellout)
    cellout{i}(1,:)=cellout{i}(1,:)./2;
    if iseven(size(x,1))
       cellout{i}(end,:)=cellout{i}(end,:)./2;
    end
end

cellout{1}=f./dt;

function[S]=avgspec(mmat1,mmat2)
eigspec=mmat1.*conj(mmat2);
S=squeeze(mean(eigspec,2));



function[s,dk]=adaptspec(eigspec,lambda,var)

s=squeeze(zeros(size(eigspec(:,1,:))));
dk=zeros(size(eigspec));
for i=1:size(eigspec,3)
    [s(:,i),dk(:,:,i)]=adaptspec_one(eigspec(:,:,i),lambda,var(i));
end
