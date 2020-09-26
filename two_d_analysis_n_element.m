% 2D ANALYSIS %%
%%program for plane stress %%with single element
E=input('Enter Youngs Modulus of the material in Pascal=');
NU=input('Enter poissons ratio=');
t=input('Enter thickness of plate in metre=');
density=input('Enter density of body in kg/m^3=');
p=input('Define plane stress(Input 1) or plane strain(Input 2) problem=');
n=input('Enter total number of node=');
K=zeros(2*n,2*n);
F=zeros(2*n,1);
e=input('Enter total number of elements=');
for l=1:e
    i=input('Enter first node of each element connectivity matrix=');
    xi=input('Enter x coordinate of first node each element in metre=');
    yi=input('Enter y coordinate of first node each element in metre=');
    j=input('Enter second node of each element connectivity matrix=');
    xj=input('Enter x coordinate of second node each element in metre=');
    yj=input('Enter y coordinate of second node each element in metre=');
    m=input('Enter third node of each element connectivity matrix=');
    xm=input('Enter x coordinate of third node each element in metre=');
    ym=input('Enter y coordinate of third node each element in metre=');
    k=LinearTriangleElementStiffness(E,NU,t,xi,yi,xj,yj,xm,ym,p);
    K=LinearTriangleAssemble(K,k,i,j,m);
    f=LinearTriangleElementBodyForceVector(density,t,xi,yi,xj,yj,xm,ym);
    F=LinearTriangleForceVectorAssemble(F,f,i,j,m);
end
pointload=zeros(2*n,1);
numberofload=input('Enter total number of point on which load will be acting=');
for i=1:numberofload
    nodenumber=input('Enter the node number on which point load will be act=');
    pointload(nodenumber,1)=input('Enter the force in that direction (consider direction also)=');
end
tractionload=zeros(2*n,1);
tractionline=input('Enter the total number of lines on which traction force will be applied=');
for i=1:tractionline
   tr1=input('Enter the first node number=');
   tr2=input('Enter the second node number=');
   t1=input('Enter the tractional load value at first node in N/m^2=');
   t2=input('Enter the tractional load value at second node in N/m^2=');
   xi=input('Enter x coordinate of first node of line=');
   yi=input('Enter y coordinate of first node of line=');
   xj=input('Enter x coordinate of first node of line=');
   yj=input('Enter y coordinate of first node of line=');
   ll=((xi-xj)^2+(yi-yj)^2)^0.5;
   c=(yj-yi)/ll;
   s=(xj-xi)/ll;
   t1x=t1*c;
   t1y=t1*s;
   t2x=t2*c;
   t2y=t2*s;
   tractionload(2*tr1-1,1)=((t*ll)/6)*(2*t1x+t2x);
   tractionload(2*tr1,1)=((t*ll)/6)*(2*t1y+t2y);
   tractionload(2*tr2-1,1)=((t*ll)/6)*(t1x+2*t2x);
   tractionload(2*tr2,1)=((t*ll)/6)*(t1y+2*t2y);
end  
Ftotal=F+pointload+tractionload;
boundary=ones(2*n,1);
m=input('Enter the number of boundary condition is known=');
for i=1:m
    q=input('Enter the number to which boundary condition is known=');
    boundary(q,1)=input('Enter the known boundary condition=');
end
disp(K);
disp(boundary);
disp(Ftotal);



