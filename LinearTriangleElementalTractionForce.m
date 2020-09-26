function [y] = LinearTriangleElementalTractionForce()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
le=((xi-xj)^2+(yi-yj)^2)^0.5;
c=(xi-xj)/(le);
s=(yi-yj)/(le);
Txi=T0i*s;
Tyi=T0i*c;
Tx=T0*s;
Ty=T0*c;

end

