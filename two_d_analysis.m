%% 2D ANALYSIS %%
%%program for plane stress %%with single element
E=input('Enter Youngs Modulus of the material in Pascal=');
NU=input('Enter poissons ratio=');
t=input('Enter thickness of plate in metre=');
p=input('Define plane stress(Input 1) or plane strain(Input 2) problem=');
n=input('Enter total number of node=');
xi=input('Enter x coordinate of first node in metre=');
yi=input('Enter y coordinate of first node in metre=');
xj=input('Enter x coordinate of second node in metre=');
yj=input('Enter y coordinate of second node in metre=');
xm=input('Enter x coordinate of third node in metre=');
ym=input('Enter y coordinate of third node in metre=');
k1=LinearTriangleElementStiffness(E,NU,t,xi,yi,xj,yj,xm,ym,p);
xi=input('Enter x coordinate of first node in metre=');
yi=input('Enter y coordinate of first node in metre=');
xj=input('Enter x coordinate of second node in metre=');
yj=input('Enter y coordinate of second node in metre=');
xm=input('Enter x coordinate of third node in metre=');
ym=input('Enter y coordinate of third node in metre=');
k2=LinearTriangleElementStiffness(E,NU,t,xi,yi,xj,yj,xm,ym,p);
K=zeros(2*n,2*n);
K=LinearTriangleAssemble(K,k1,1,3,4);
K=LinearTriangleAssemble(K,k1,1,2,3);
disp(K);