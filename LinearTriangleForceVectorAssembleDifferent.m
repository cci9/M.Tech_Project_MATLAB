function [y] = LinearTriangleForceVectorAssembleDifferent(F,f,i,j,m,PointLoad)
F(2*i-1,1) = F(2*i-1,1) + f(1,1)+PointLoad(2*i-1,1);
F(2*i,1) = F(2*i,1) + f(2,1)+PointLoad(2*i,1);
F(2*j-1,1) = F(2*j-1,1) + f(3,1)+PointLoad(2*i-1,1);
F(2*j,1) = F(2*j,1) + f(4,1)+PointLoad(2*i,1);
F(2*m-1,1) = F(2*m-1,1) + f(5,1)+PointLoad(2*i-1,1);
F(2*m,1) = F(2*m,1) + f(6,1)+PointLoad(2*i,1);
y=F;
end

