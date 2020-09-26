function [y] = LinearTriangleElementArea(xi,yi,xj,yj,xm,ym)
y=(xi*(yj-ym)+xj*(ym-yi)+xm*(yi-yj))/2;
end

