% 2D ANALYSIS %%
%%program for plane stress %%with single element
E=input('Enter Youngs Modulus of the material in Pascal=');
NU=input('Enter poissons ratio=');
t=input('Enter thickness of plate in metre=');
Fx=input('Enter x direction body force (+ve for right directional load and -ve for left directional load) in Newton=');
Fy=input('Enter y direction body force (+ve for upward directional load and -ve for downward directional load) in Newton=');
p=input('Define plane stress(Input 1) or plane strain(Input 2) problem=');
n=input('Enter total number of node=');
K=zeros(2*n,2*n);
F=zeros(2*n,1);
e=input('Enter total number of elements=');
for l=1:e
    i=input('Enter first node of each element connectivity matrix=');
    xi=input('Enter x coordinate of first node each element in metre=');
    yi=input('Enter y coordinate of first node each element in metre=');
    Pix=input('Enter the point load acting at first node in x direction in newton=');
    Piy=input('Enter the point load acting at first node in y direction in newton=');
    j=input('Enter second node of each element connectivity matrix=');
    xj=input('Enter x coordinate of second node each element in metre=');
    yj=input('Enter y coordinate of second node each element in metre=');
    Pjx=input('Enter the point load acting at second node in x direction in newton=');
    Pjy=input('Enter the point load acting at second node in y direction in newton=');
    m=input('Enter third node of each element connectivity matrix=');
    xm=input('Enter x coordinate of third node each element in metre=');
    ym=input('Enter y coordinate of third node each element in metre=');
    Pmx=input('Enter the point load acting at third node in x direction in newton=');
    Pmy=input('Enter the point load acting at third node in y direction in newton=');
    k=LinearTriangleElementStiffness(E,NU,t,xi,yi,xj,yj,xm,ym,p);
    K=LinearTriangleAssemble(K,k,i,j,m);
    f=LinearTriangleElementBodyForceVector(Fx,Fy,t,xi,yi,xj,yj,xm,ym);
    F=LinearTriangleForceVectorAssemble(F,f,i,j,m,PointLoad);
end
for o=1:n
    i=input('Enter the node number at which load is to be added in newton=');
    Px=input('Enter the point load acting at node in x direction in newton=');
    Py=input('Enter the point load acting at node in y direction in newton=');
    PointLoad=LinearTriangleElementalConcentratedLoad(Pix,Piy,Pjx,Pjy,Pmx,Pmy);
    
disp(K);
disp(F);
