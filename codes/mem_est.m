function [mem, x, y] = mem_est(a1,a2,b1,b2)
%% Bill's code

d1=a1;d2=b1;d3=a2;d4=b2;

c1=complex(1.,0).*d1+complex(0,1.).*d2;
c2=complex(1.,0).*d3+complex(0,1.).*d4;

p1=(c1-c2.*conj(c1))./(1-abs(c1).^2);
p2=c2-c1.*p1;

x=1.-p1.*conj(c1)-p2.*conj(c2);

% calculate MEM using 1 deg resolution
a=1:360;
a=a.*pi/180;
e1=complex(1.,0)*cos(a)-complex(0,1.)*sin(a);
e2=complex(1.,0)*cos(2*a)-complex(0,1.)*sin(2*a);
y=abs(complex(1.,0)-p1*e1-p2*e2).^2;

mem=abs((x*ones(1,360))./y);

% normalize MEM distribution
mem=((1./sum(mem'))'*ones(1,360)).*mem;

%  direction indices are reveresed after matrix
%  calcs, so flip directions back to native 
%  direction coordinates.
mem=fliplr(mem);
return
