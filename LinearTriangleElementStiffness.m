function [y] = LinearTriangleElementStiffness(E,NU,t,xi,yi,xj,yj,xm,ym,p)
Ae=(xi*(yj-ym)+xj*(ym-yi)+xm*(yi-yj))/2;
b11i=yj-ym;
b11j=ym-yi;
b11m=yi-yj;
b12i=xm-xj;
b12j=xi-xm;
b12m=xj-xi;
B=(1/(2*Ae))*[b11i 0 b11j 0 b11m 0;0 b12i 0 b12j 0 b12m;b12i b11i b12j b11j b12m b11m];
if p==1
    En=E;
    NUn=NU;
    D=(En/(1-NUn*NUn))*[1 NUn 0;NUn 1 0;0 0 (1-NUn)/2];
elseif p==2
    En=E/(1-NU*NU);
    NUn=NU/(1-NU);
    D=(En/(1-NUn*NUn))*[1 NUn 0;NUn 1 0;0 0 (1-NUn)/2];
end
y=t*Ae*B'*D*B;

