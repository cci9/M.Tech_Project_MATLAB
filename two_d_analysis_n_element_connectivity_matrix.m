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
K=zeros(2*n,2*n);
for i=1:n
    xy(i,1)=input('Enter the x coordinate of the each node in sequence =');
    xy(i,2)=input('Enter the y coordinate of the each node in sequence =');
end
for i=1:e
    xn=zeros(1,3);
    yn=zeros(1,3);
    for j=1:3
        xn(1,j)=xy(cm(i,j),1);
        yn(1,j)=xy(cm(i,j),2);    
    end
    Ae(i)=(xn(1,1)*(yn(1,2)-yn(1,3))+xn(1,2)*(yn(1,3)-yn(1,1))+xn(1,3)*(yn(1,1)-yn(1,2)))*0.5;
    Ae(i)=abs(Ae(i));
    B=(0.5/(Ae(i)))*[yn(1,2)-yn(1,3) 0 yn(1,3)-yn(1,1) 0 yn(1,1)-yn(1,2) 0;
        0 xn(1,3)-xn(1,2) 0 xn(1,1)-xn(1,3) 0 xn(1,2)-xn(1,1);
        xn(1,3)-xn(1,2) yn(1,2)-yn(1,3) xn(1,1)-xn(1,3) yn(1,3)-yn(1,1) xn(1,2)-xn(1,1) yn(1,1)-yn(1,2)]; 
    k=t*Ae(i)*B'*D*B;
    cm1=cm(i,1);
    cm2=cm(i,2);
    cm3=cm(i,3);  
    K([2*cm1-1:2*cm1,2*cm2-1:2*cm2,2*cm3-1:2*cm3],[2*cm1-1:2*cm1,2*cm2-1:2*cm2,2*cm3-1:2*cm3])=K([2*cm1-1:2*cm1,2*cm2-1:2*cm2,2*cm3-1:2*cm3],[2*cm1-1:2*cm1,2*cm2-1:2*cm2,2*cm3-1:2*cm3])+k;
end
disp('GLOBAL STIFFNESS MATRIX');
disp(K);
fbn=zeros(2*n,1);
fb=zeros(2*n,1);
for i=1:e
     for j=1:3
        ii=cm(i,j);
        fbn(2*ii-1,1)=0;
        fbn(2*ii,1)=-(density*9.81/3)*Ae(i)*t;
     end
     fb=fb+fbn;
end
ft=zeros(2*n,1);
tractionline=input('Enter the total number of lines on which traction force will be applied=');
for i=1:tractionline
   tr1=input('Enter the first node number on traction force will act=');
   tr2=input('Enter the second node number on traction force will act=');
   t1=input('Enter the tractional load value at first node in N/m^2=');
   t2=input('Enter the tractional load value at second node in N/m^2=');
   x1=xy(tr1,1);
   y1=xy(tr1,2);
   x2=xy(tr2,1);
   y2=xy(tr2,2);
   ll=sqrt((x2-x1)^2+(y2-y1)^2); % length of elemental line on which tractional load wil act
   s=(y1-y2)/ll; %angle of line on which tractional load will be applied making angle with horizontal
   c=(x2-x1)/ll; %angle of line on which tractional load will be applied making angle with horizontal
   tx1=s*t1;     %resolved component of tractional load on first node
   ty1=c*t1;     %resolved component of tractional load on first node
   tx2=s*t2;     %resolved component of tractional load on second node
   ty2=c*t2;     %resolved component of tractional load on second node
   ft(2*tr1-1,1)=(t*ll/6)*(2*tx1+tx2);
   ft(2*tr1,1)=(t*ll/6)*(2*ty1+ty2);
   ft(2*tr2-1,1)=(t*ll/6)*(tx1+2*tx2);
   ft(2*tr2,1)=(t*ll/6)*(ty1+2*ty2);
end
fp=zeros(2*n,1);
numberofload=input('Enter total number of point on which load will be acting=');
for i=1:numberofload
    nodenumber=input('Enter the row number on which point load will be act(2*node-1 for Fx component and 2*node for Fy component)=');
    fp(nodenumber,1)=input('Enter the force with considering +ve or -ve direction in Newton=');
end
F=fb+ft+fp;
disp('GLOBAL LOAD VECTOR');
disp(F);
% Lets apply the boundary condition
u=ones(2*n,1);
uknown=input('Enter total number of boundary condition are known=');
for i=1:uknown
    uknownpoint=input('Enter the row number of which boundary condition is known(2*node-1 for ux disp. and 2*node for uy disp.)=');
    u(uknownpoint,1)=input('Enter the boundary condition with (considering +ve or -ve direction) in metre=');
end
% Lets reduce the matrix size
Kred=K;
for i=1:2*n
    if(u(i)==0)
      if(i==1)
      Kred=Kred((2:2*n),(2:2*n));
      end
      if(i>1)
          j=i-(length(K)-length(Kred));
          [m,n] = size(Kred);
          Kred1 = Kred(1:j-1,1:j-1);
          Kred2 = Kred(1:j-1,j+1:n);
          Kred3 = Kred(j+1:m,1:j-1);
          Kred4 = Kred(j+1:m,j+1:n);
          Kred = [Kred1 Kred2; Kred3 Kred4];
      end
    end
end
disp('REDUCED STIFFNESS VECTOR AFTER APPLYING BOUNDARY CONDITION');
disp(Kred);
[a,b]=size(F);
Fred=F;
[c,d]=size(Fred);
n=input('Enter total number of nodes again=');
for i=1:2*n
    if(u(i)==0)
        if(i==1)
            Fred=F((2:2*n),1);
            [c,d]=size(Fred);
        end
        if(i>1)
            [c,d]=size(Fred);
            j=i-(a-c);
            Fred1=Fred(1:j-1,d);
            Fred2=Fred(j+1:c,d);
            Fred=[Fred1;Fred2];
        end
    end
end
disp('REDUCED FORCE VECTOR AFTER APPLYING BOUNDARY CONDITION');
disp(Fred);
unew=(inv(Kred))*Fred;
disp('DEGREE OF FREEDOM FOR AFTER APPLYING BOUNDARY CONDION');
disp(unew);
% Lets find Reaction on the fixed supports
j=1;
for i=1:2*n
        if (u(i)~=0)
            u(i)=unew(j);
            j=j+1;
        end
end
disp('DEGREE OF FREEDOM FOR ALL NODE');
disp(u);
Reactions=(K*u)-F;
disp('REACTIONS IN VECTOR FORM');
disp(Reactions);
ustress=zeros(6,1);
for i=1:e
    xn=zeros(1,3);
    yn=zeros(1,3);
    for j=1:3
        xn(1,j)=xy(cm(i,j),1);
        yn(1,j)=xy(cm(i,j),2);    
    end
    Ae=(xn(1,1)*(yn(1,2)-yn(1,3))+xn(1,2)*(yn(1,3)-yn(1,1))+xn(1,3)*(yn(1,1)-yn(1,2)))*0.5;
    Ae=abs(Ae);
    B=(0.5/Ae)*[yn(1,2)-yn(1,3) 0 yn(1,3)-yn(1,1) 0 yn(1,1)-yn(1,2) 0;
        0 xn(1,3)-xn(1,2) 0 xn(1,1)-xn(1,3) 0 xn(1,2)-xn(1,1);
        xn(1,3)-xn(1,2) yn(1,2)-yn(1,3) xn(1,1)-xn(1,3) yn(1,3)-yn(1,1) xn(1,2)-xn(1,1) yn(1,1)-yn(1,2)];
    
    m=1;
    for j=1:3
        nd=input('Enter the nodes of each connectivity matrix for getting the stresses in anticlockwise sequence=');
        ustress(m,1)=u(2*nd-1);
        ustress(m+1,1)=u(2*nd);
        m=m+2;
    end
    disp('DEGREE OF FREEDOM FOR EACH ELEMENT');
    disp(ustress);
    StrainElement=B*ustress;
    StressElement=D*B*ustress;
    disp('STRAIN IN EACH ELEMENT');
    disp(StrainElement);
    disp('STRESSS IN EACH ELEMENT');
    disp(StressElement);
    sigma1=0.5*(StressElement(1)+StressElement(2))+((0.5*(StressElement(1)-StressElement(2)))^2+(StressElement(3))^2)^0.5;
    sigma2=0.5*(StressElement(1)+StressElement(2))-((0.5*(StressElement(1)-StressElement(2)))^2+(StressElement(3))^2)^0.5;
    sigma_vonmis(i)=(sigma1^2+sigma2^2-sigma1*sigma2)^0.5;
    disp('VONMISES STRESSES FOR EACH ELEMENT');
    disp(sigma_vonmis(i));
end
disp('Element number\t Elemental Von-Mises Stress');
for i=1:e
    fprintf('%2.0f \t %8.3f \n',i,sigma_vonmis(i));
end

