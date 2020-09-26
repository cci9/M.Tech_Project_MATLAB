% 2d analysis by  quadratic trianguar element(LST) %
E=210*10^9;% E=input('Enter Youngs Modulus of the material in Pascal=');      
NU=0.3;% NU=input('Enter poissons ratio=');
p=1;% p=input('Define plane stress problem(Input 1) or define plane strain problem(Input 2)=');
t=0.025;% t=input('Enter thickness of plate in metre=');
density=0;% density=input('Enter density of body in kg/m^3=');
n=6;% n=input('Enter total number of node=');
e=1;% e=input('Enter total number of elements=');
if p==1
    En=E;
    NUn=NU;
    D=(En/(1-NUn*NUn))*[1 NUn 0;NUn 1 0;0 0 (1-NUn)/2];
elseif p==2
    En=E/(1-NU*NU);
    NUn=NU/(1-NU);
    D=(En/(1-NUn*NUn))*[1 NUn 0;NUn 1 0;0 0 (1-NUn)/2];
end
cm=zeros(e,6);
for i=1:e
    for j=1:6
        cm(i,j)=input('Enter the nodes of each connectivity matrix in anticlockwise sequence=');
    end
end
for i=1:n
    xy(i,1)=input('Enter the x coordinate of the each node in sequence =');
    xy(i,2)=input('Enter the y coordinate of the each node in sequence =');
end
K=zeros(2*n,2*n);
Kelemental=zeros(12,12);
for i=1:e
    xn=zeros(1,6);
    yn=zeros(1,6);
    for j=1:6
        xn(1,j)=xy(cm(i,j),1);
        yn(1,j)=xy(cm(i,j),2);    
    end
    for y=1:7
        if y==1
           weight=3/120;
           g1=1;
           g2=0;
        end
        if y==2
           weight=3/120;
           g1=0;
           g2=1;
        end
        if y==3
           weight=3/120;
           g1=0;
           g2=0;
        end
        if y==4
           weight=8/120;
           g1=0.5;
           g2=0.5;
        end
        if y==5
           weight=8/120;
           g1=0.5;
           g2=0;
        end
        if y==6
           weight=8/120;
           g1=0;
           g2=0.5;
        end
        if y==7
           weight=27/120;
           g1=1/3;
           g2=1/3;
        end
        do=[4*g1-1 0 4*g1+4*g2-3 4*g2 -4*g2 4-8*g1-4*g2;0 4*g1-1 4*g1+4*g2-3 4*g1 4-4*g1-8*g2 -4*g1];
        xylocal=[xn(1,1) yn(1,1);xn(1,2) yn(1,2);xn(1,3) yn(1,3);xn(1,4) yn(1,4);xn(1,5) yn(1,5);xn(1,6) yn(1,6)];
        j=do*xylocal;
        G=(inv(j))*do;
        B=[G(1,1) 0 G(1,2) 0 G(1,3) 0 G(1,4) 0 G(1,5) 0 G(1,6) 0;0 G(2,1) 0 G(2,2) 0 G(2,3) 0 G(2,4) 0 G(2,5) 0 G(2,6);G(2,1) G(1,1) G(2,2) G(1,2) G(2,3) G(1,3) G(2,4) G(1,4) G(2,5) G(1,5) G(2,6) G(1,6)];
        Kgauss=(B')*D*B*det(j)*t*weight;
        Kelemental=Kelemental+Kgauss;
    end
    cm1=cm(i,1);
    cm2=cm(i,2);
    cm3=cm(i,3); 
    cm4=cm(i,4);
    cm5=cm(i,5);
    cm6=cm(i,6);
    K([2*cm1-1:2*cm1,2*cm2-1:2*cm2,2*cm3-1:2*cm3,2*cm4-1:2*cm4,2*cm5-1:2*cm5,2*cm6-1:2*cm6],[2*cm1-1:2*cm1,2*cm2-1:2*cm2,2*cm3-1:2*cm3,2*cm4-1:2*cm4,2*cm5-1:2*cm5,2*cm6-1:2*cm6])=K([2*cm1-1:2*cm1,2*cm2-1:2*cm2,2*cm3-1:2*cm3,2*cm4-1:2*cm4,2*cm5-1:2*cm5,2*cm6-1:2*cm6],[2*cm1-1:2*cm1,2*cm2-1:2*cm2,2*cm3-1:2*cm3,2*cm4-1:2*cm4,2*cm5-1:2*cm5,2*cm6-1:2*cm6])+Kelemental;
end
disp('GLOBAL STIFFNESS MATRIX');
disp(K);
fb=zeros(2*n,1);
fbn=zeros(2*n,1);
for i=1:e
     for j=1:6
        ii=cm(i,j);
        fb(2*ii-1,1)=0;
        if j==1
            N=g1*(2*g1-1);
        end
        if j==2
            N=g2*(2*g2-1);
        end
        if j==3
            N=(1-g1-g2)*(2*(1-g1-g2)-1);
        end
        if j==4
            N=4*g1*g2;
        end
        if j==5
            N=4*g2*(1-g1-g2);
        end
        if j==6
            N=4*g1*(1-g1-g2);
        end
        for x=1:7
            if x==1
               weight=3/120;
               g1=1;
               g2=0;
            end
            if x==2
               weight=3/120;
               g1=0;
               g2=1;
            end
            if x==3
               weight=3/120;
               g1=0;
               g2=0;
            end
            if x==4
               weight=8/120;
               g1=0.5;
               g2=0.5;
            end
            if x==5
               weight=8/120;
               g1=0.5;
               g2=0;
            end
            if x==6
               weight=8/120;
               g1=0;
               g2=0.5;
            end
            if x==7
               weight=27/120;
               g1=1/3;
               g2=1/3;
            end   
            do=[4*g1-1 0 4*g1+4*g2-3 4*g2 -4*g2 4-8*g1-4*g2;0 4*g1-1 4*g1+4*g2-3 4*g1 4-4*g1-8*g2 -4*g1];
            xylocal=[xn(1,1) yn(1,1);xn(1,2) yn(1,2);xn(1,3) yn(1,3);xn(1,4) yn(1,4);xn(1,5) yn(1,5);xn(1,6) yn(1,6)];
            j=do*xylocal;
            fb(2*ii,1)=-(density*9.81)*N*det(j)*t;
        end
     end
     fbn=fbn+fb;
end
ft=zeros(2*n,1);
tractionline=input('Enter the total number of lines on which traction force will be applied=');
for i=1:tractionline
   tr1=input('Enter the first node number in anticlockwise sequence on traction force will act=');
   tr2=input('Enter the second node number in anticlockwise sequence on traction force will act=');
   tr3=input('Enter the third node number in anticlockwise sequence on traction force will act=');
   T=input('Enter the UDL tractional load value on line in N/m^2=');
   x1=xy(tr1,1);
   y1=xy(tr1,2);
   x2=xy(tr3,1);
   y2=xy(tr3,2);
   ll=sqrt((x2-x1)^2+(y2-y1)^2); % length of elemental line on which tractional load wil act
   s=(y1-y2)/ll; %angle of line on which tractional load will be applied making angle with horizontal
   c=(x2-x1)/ll; %angle of line on which tractional load will be applied making angle with horizontal
   ft(2*tr1-1,1)=(t*ll/6)*s*T;
   ft(2*tr1,1)=(t*ll/6)*c*T;
   ft(2*tr2-1,1)=(2*t*ll/3)*s*T;
   ft(2*tr2,1)=(2*t*ll/3)*c*T;
   ft(2*tr3-1,1)=(t*ll/6)*s*T;
   ft(2*tr3,1)=(t*ll/6)*c*T;
end
fp=zeros(2*n,1);
numberofload=input('Enter total number of point on which load will be acting=');
for i=1:numberofload
    nodenumber=input('Enter the row number on which point load will be act(2*node-1 for Fx component and 2*node for Fy component)=');
    fp(nodenumber,1)=input('Enter the force with considering +ve or -ve direction in Newton=');
end
F=fbn+ft+fp;
disp('GLOBAL LOAD VECTOR');
disp(F);
% % Lets apply the boundary condition
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
    if u(i)==0
        if i==1
            Fred=F((2:2*n),1);
            [c,d]=size(Fred);
        end
        if i>1
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
    xn=zeros(1,6);
    yn=zeros(1,6);
    for j=1:6
        xn(1,j)=xy(cm(i,j),1);
        yn(1,j)=xy(cm(i,j),2);    
    end
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
    end
    do=[4*g1-1 0 4*g1+4*g2-3 4*g2 -4*g2 4-8*g1-4*g2;0 4*g1-1 4*g1+4*g2-3 4*g1 4-4*g1-8*g2 -4*g1];
    xylocal=[xn(1,1) yn(1,1);xn(1,2) yn(1,2);xn(1,3) yn(1,3);xn(1,4) yn(1,4);xn(1,5) yn(1,5);xn(1,6) yn(1,6)];
    j=do*xylocal;
    G=(inv(j))*do;
    B=[G(1,1) 0 G(1,2) 0 G(1,3) 0 G(1,4) 0 G(1,5) 0 G(1,6) 0;0 G(2,1) 0 G(2,2) 0 G(2,3) 0 G(2,4) 0 G(2,5) 0 G(2,6);G(2,1) G(1,1) G(2,2) G(1,2) G(2,3) G(1,3) G(2,4) G(1,4) G(2,5) G(1,5) G(2,6) G(1,6)];
    m=1;
    for j=1:6
        ustress(m,1)=u(2*cm(i,j)-1);
        ustress(m+1,1)=u(2*cm(i,j));
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
disp('Element number      Elemental Von-Mises Stress');
for i=1:e
    fprintf('%2.0f \t\t\t %8.3f \n',i,sigma_vonmis(i));
end
vms=max(sigma_vonmis);
disp(vms);
for i=1:e
    xn=zeros(1,3);
    yn=zeros(1,3);
    for j=1:3
        xn(1,j)=xy(cm(i,j),1);
        yn(1,j)=xy(cm(i,j),2);    
    end
    x1=[xn(1,1),xn(1,2)];
    x2=[xn(1,2),xn(1,3)];
    x3=[xn(1,3),xn(1,1)];
    y1=[yn(1,1),yn(1,2)];
    y2=[yn(1,2),yn(1,3)];
    y3=[yn(1,3),yn(1,1)];
    line(x1,y1);
    hold on
    line(x2,y2);
    hold on
    line(x3,y3);
    hold on
    if (sigma_vonmis(i))>0.95*vms
        fill([xn(1,1) xn(1,2) xn(1,3)],[yn(1,1) yn(1,2) yn(1,3)],'r');
        hold on
    end
    if 0.8*vms<(sigma_vonmis(i))<0.95*vms
        fill([xn(1,1) xn(1,2) xn(1,3)],[yn(1,1) yn(1,2) yn(1,3)],'y');
        hold on
    end
    if 0.6*vms<(sigma_vonmis(i))<0.8*vms
        fill([xn(1,1) xn(1,2) xn(1,3)],[yn(1,1) yn(1,2) yn(1,3)],'g');
        hold on
    end
    if 0.4*vms<(sigma_vonmis(i))<0.6*vms
        fill([xn(1,1) xn(1,2) xn(1,3)],[yn(1,1) yn(1,2) yn(1,3)],'c');
        hold on
    end
    if (sigma_vonmis(i))<0.4*vms
        fill([xn(1,1) xn(1,2) xn(1,3)],[yn(1,1) yn(1,2) yn(1,3)],'b');
        hold on
    end
end 
