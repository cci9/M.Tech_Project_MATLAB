% 2d analysis by  trianguar element %
E=210*10^9;% E=input('Enter Youngs Modulus of the material in Pascal=');      
NU=0.3;% NU=input('Enter poissons ratio=');
p=1;% p=input('Define plane stress problem(Input 1) or define plane strain problem(Input 2)=');

t=0.025;% t=input('Enter thickness of plate in metre=');
density=2000;% density=input('Enter density of body in kg/m^3=');
n=4;% n=input('Enter total number of node=');
e=2;% e=input('Enter total number of elements=');
if p==1
    En=E;
    NUn=NU;
    D=(En/(1-NUn*NUn))*[1 NUn 0;NUn 1 0;0 0 (1-NUn)/2];
elseif p==2
    En=E/(1-NU*NU);
    NUn=NU/(1-NU);
    D=(En/(1-NUn*NUn))*[1 NUn 0;NUn 1 0;0 0 (1-NUn)/2];
end
cm=zeros(e,3);
for i=1:e
    for j=1:3
        cm(i,j)=input('Enter the nodes of each connectivity matrix in anticlockwise sequence=');
    end
end
for i=1:n
    xy(i,1)=input('Enter the x coordinate of the each node in sequence =');
    xy(i,2)=input('Enter the y coordinate of the each node in sequence =');
end
fb=zeros(2*n,1);
s=zeros(2*n,1);
K=zeros(2*n,2*n);
Kn=zeros(2*n,2*n);
for i=1:e
    xn=zeros(1,3);
    yn=zeros(1,3);
    for j=1:3
        xn(1,j)=xy(cm(i,j),1);
        yn(1,j)=xy(cm(i,j),2);    
    end
    Ae(i)=(xn(1,1)*(yn(1,2)-yn(1,3))+xn(1,2)*(yn(1,3)-yn(1,1))+xn(1,3)*(yn(1,1)-yn(1,2)))*0.5;
    Ae(i)=abs(Ae(i));
    disp(Ae(i));
    B=(0.5/(Ae(i)))*[yn(1,2)-yn(1,3) 0 yn(1,3)-yn(1,1) 0 yn(1,1)-yn(1,2) 0;
        0 xn(1,3)-xn(1,2) 0 xn(1,1)-xn(1,3) 0 xn(1,2)-xn(1,1);
        xn(1,3)-xn(1,2) yn(1,2)-yn(1,3) xn(1,1)-xn(1,3) yn(1,3)-yn(1,1) xn(1,2)-xn(1,1) yn(1,1)-yn(1,2)]; 
    k=t*Ae(i)*B'*D*B;
    cm1=cm(i,1);
    cm2=cm(i,2);
    cm3=cm(i,3);  
    K(2*cm1-1,2*cm1-1) = K(2*cm1-1,2*cm1-1) + k(1,1);
    K(2*cm1-1,2*cm1) = K(2*cm1-1,2*cm1) + k(1,2);
    K(2*cm1-1,2*cm2-1) = K(2*cm1-1,2*cm2-1) + k(1,3);
    K(2*cm1-1,2*cm2) = K(2*cm1-1,2*cm2) + k(1,4);
    K(2*cm1-1,2*cm3-1) = K(2*cm1-1,2*cm3-1) + k(1,5);
    K(2*cm1-1,2*cm3) = K(2*cm1-1,2*cm3) + k(1,6);
    K(2*cm1,2*cm1-1) = K(2*cm1,2*cm1-1) + k(2,1);
    K(2*cm1,2*cm1) = K(2*cm1,2*cm1) + k(2,2);
    K(2*cm1,2*cm2-1) = K(2*cm1,2*cm2-1) + k(2,3);
    K(2*cm1,2*cm2) = K(2*cm1,2*cm2) + k(2,4);
    K(2*cm1,2*cm3-1) = K(2*cm1,2*cm3-1) + k(2,5);
    K(2*cm1,2*cm3) = K(2*cm1,2*cm3) + k(2,6);
    K(2*cm2-1,2*cm1-1) = K(2*cm2-1,2*cm1-1) + k(3,1);
    K(2*cm2-1,2*cm1) = K(2*cm2-1,2*cm1) + k(3,2);
    K(2*cm2-1,2*cm2-1) = K(2*cm2-1,2*cm2-1) + k(3,3);
    K(2*cm2-1,2*cm2) = K(2*cm2-1,2*cm2) + k(3,4);
    K(2*cm2-1,2*cm3-1) = K(2*cm2-1,2*cm3-1) + k(3,5);
    K(2*cm2-1,2*cm3) = K(2*cm2-1,2*cm3) + k(3,6);
    K(2*cm2,2*cm1-1) = K(2*cm2,2*cm1-1) + k(4,1);
    K(2*cm2,2*cm1) = K(2*cm2,2*cm1) + k(4,2);
    K(2*cm2,2*cm2-1) = K(2*cm2,2*cm2-1) + k(4,3);
    K(2*cm2,2*cm2) = K(2*cm2,2*cm2) + k(4,4);
    K(2*cm2,2*cm3-1) = K(2*cm2,2*cm3-1) + k(4,5);
    K(2*cm2,2*cm3) = K(2*cm2,2*cm3) + k(4,6);
    K(2*cm3-1,2*cm1-1) = K(2*cm3-1,2*cm1-1) + k(5,1);
    K(2*cm3-1,2*cm1) = K(2*cm3-1,2*cm1) + k(5,2);
    K(2*cm3-1,2*cm2-1) = K(2*cm3-1,2*cm2-1) + k(5,3);
    K(2*cm3-1,2*cm2) = K(2*cm3-1,2*cm2) + k(5,4);
    K(2*cm3-1,2*cm3-1) = K(2*cm3-1,2*cm3-1) + k(5,5);
    K(2*cm3-1,2*cm3) = K(2*cm3-1,2*cm3) + k(5,6);
    K(2*cm3,2*cm1-1) = K(2*cm3,2*cm1-1) + k(6,1);
    K(2*cm3,2*cm1) = K(2*cm3,2*cm1) + k(6,2);
    K(2*cm3,2*cm2-1) = K(2*cm3,2*cm2-1) + k(6,3);
    K(2*cm3,2*cm2) = K(2*cm3,2*cm2) + k(6,4);
    K(2*cm3,2*cm3-1) = K(2*cm3,2*cm3-1) + k(6,5);
    K(2*cm3,2*cm3) = K(2*cm3,2*cm3) + k(6,6);
%     Kn(2*e-1:2*e,2*e-1:2*e)=Kn(2*e-1:2*e,2*e-1:2*e)+k;
    Kn([2*cm1-1:2*cm1,2*cm2-1:2*cm2,2*cm3-1:2*cm3],[2*cm1-1:2*cm1,2*cm2-1:2*cm2,2*cm3-1:2*cm3])=Kn([2*cm1-1:2*cm1,2*cm2-1:2*cm2,2*cm3-1:2*cm3],[2*cm1-1:2*cm1,2*cm2-1:2*cm2,2*cm3-1:2*cm3])+k;
end
disp(K);
 disp(Kn);
% for i=1:e
%     for j=1:3
%         ii=cm(i,j);
%         Kn(2*ii-1:2*ii,2*
%     end
% end