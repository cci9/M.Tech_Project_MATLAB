function [y] = LinearTriangleElementBodyForceVector(density,t,xi,yi,xj,yj,xm,ym)
Ae=(xi*(yj-ym)+xj*(ym-yi)+xm*(yi-yj))/2;
y=(Ae*t/3)*[0;-(density*9.81);0;-(density*9.81);0;-(density*9.81)];
end

