%% 2D ANALYSIS %%
%%program for plane stress %%with single element
E=input('Enter Youngs Modulus of the material in Pascal=');
mu=input('Enter poissons ratio=');
t=input('Enter thickness of plate in metre=');
for i=1:3
    x(i)=input('Enter x coordinate of node in metre=');
    y(i)=input('Enter y coordinate of node in metre=');
end
Ae=0.5*(x(1)*(y(2)-y(3))+x(2)*(y(3)-y(1))+x(3)*(y(1)-y(2)));
B=(1/2*Ae)*[y(2)-y(3) 0 y(3)-y(1) 0 y(1)-y(2) 0;0 x(3)-x(2) 0 x(1)-x(3) 0 x(2)-x(1);x(3)-x(2) y(2)-y(3) x(1)-x(3) y(3)-y(1) x(2)-x(1) y(1)-y(2)];
D=(E/(1-mu*mu))*[1 mu 0;mu 1 0;0 0 (1-mu)/2];
disp(Ae);
disp(D);
disp(B);
disp(B');
K=(B'*D*B)*Ae*t;
disp(K);

