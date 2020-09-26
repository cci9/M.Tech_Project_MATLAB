%%two d analysis linear starin triangle one element
E=210*10^9;% E=input('Enter Youngs Modulus of the material in Pascal=');      
NU=0.3;% NU=input('Enter poissons ratio=');
p=1;% p=input('Define plane stress problem(Input 1) or define plane strain problem(Input 2)=');
t=1;% t=input('Enter thickness of plate in metre=');
density=0;% density=input('Enter density of body in kg/m^3=');
% n=4;% n=input('Enter total number of node=');
% e=2;% e=input('Enter total number of elements=');
x1=4;
y1=1;
x2=1;
y2=4;
x3=1;
y3=1;
x4=2.5;
y4=2.5;
x5=1;
y5=2.5;
x6=2.5;
y6=1;
if p==1
    En=E;
    NUn=NU;
    D=(En/(1-NUn*NUn))*[1 NUn 0;NUn 1 0;0 0 (1-NUn)/2];
elseif p==2
    En=E/(1-NU*NU);
    NUn=NU/(1-NU);
    D=(En/(1-NUn*NUn))*[1 NUn 0;NUn 1 0;0 0 (1-NUn)/2];
end
Kelemental=zeros(12,12);
for n=1:7
    if n==1
        weight=3/120;
        g1=1;
        g2=0;
    end
    if n==2
        weight=3/120;
        g1=0;
        g2=1;
    end
    if n==3
        weight=3/120;
        g1=0;
        g2=0;
    end
    if n==4
        weight=8/120;
        g1=0.5;
        g2=0.5;
    end
    if n==5
        weight=8/120;
        g1=0.5;
        g2=0;
    end
    if n==6
        weight=8/120;
        g1=0;
        g2=0.5;
    end
    if n==7
        weight=27/120;
        g1=1/3;
        g2=1/3;
    end
do=[4*g1-1 0 4*g1+4*g2-3 4*g2 -4*g2 4-8*g1-4*g2;0 4*g1-1 4*g1+4*g2-3 4*g1 4-4*g1-8*g2 -4*g1];
xylocal=[x1 y1;x2 y2;x3 y3;x4 y4;x5 y5;x6 y6];
j=do*xylocal;
G=(inv(j))*do;
disp(G);
B=[G(1,1) 0 G(1,2) 0 G(1,3) 0 G(1,4) 0 G(1,5) 0 G(1,6) 0;0 G(2,1) 0 G(2,2) 0 G(2,3) 0 G(2,4) 0 G(2,5) 0 G(2,6);G(2,1) G(1,1) G(2,2) G(1,2) G(2,3) G(1,3) G(2,4) G(1,4) G(2,5) G(1,5) G(2,6) G(1,6)];
Kgauss=(B')*D*B*det(j)*t*weight;
Kelemental=Kelemental+Kgauss;
% disp(Kelemental);
end
disp(Kelemental);


