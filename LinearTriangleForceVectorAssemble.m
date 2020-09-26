function [y] = LinearTriangleForceVectorAssemble(F,f,i,j,m)
F(2*i-1,1) = F(2*i-1,1) + f(1,1);
F(2*i,1) = F(2*i,1) + f(2,1);
F(2*j-1,1) = F(2*j-1,1) + f(3,1);
F(2*j,1) = F(2*j,1) + f(4,1);
F(2*m-1,1) = F(2*m-1,1) + f(5,1);
F(2*m,1) = F(2*m,1) + f(6,1);
y=F;
end

