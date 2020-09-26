 %%%%%%%%% Path Generation Four Precision Point Problem without Prescribed Timing %%%%%%%%%%
% Fixed line A0B0
a0x=-2;%a0x=input('Enter x coordinate of A0 point=');
a0y=2;%a0y=input('Enter y coordinate of A0 point=');

% Different Position
p1x=3;%p1x=input('Enter x coordinate of P1 point=');
p1y=8;%p1y=input('Enter y coordinate of P1 point=');

p2x=4;%p2x=input('Enter x coordinate of P2 point=');
p2y=8.3;%p2y=input('Enter y coordinate of P2 point=');

p3x=5;%p3x=input('Enter x coordinate of P3 point=');
p3y=8.5;%p3y=input('Enter y coordinate of P3 point=');

p4x=6;%p3x=input('Enter x coordinate of P3 point=');
p4y=8.3;%p3y=input('Enter y coordinate of P3 point=');

plot(p1x,p1y,'r*','MarkerSize',3);
hold on
plot(p2x,p2y,'r*','MarkerSize',3);
hold on
plot(p3x,p3y,'r*','MarkerSize',3);
hold on
plot(p4x,p4y,'r*','MarkerSize',3);
hold on
plot(a0x,a0y,'r*','MarkerSize',3);
hold on

% Finding perpendicular bisector of P1P4
line([p1x p4x],[p1y p4y]);
hold on
p1p4xCenter=(p1x+p4x)/2;
p1p4yCenter=(p1y+p4y)/2;
plot(p1p4xCenter, p1p4yCenter, 'r*', 'MarkerSize', 3);
hold on
mp1p4=(p1y-p4y)/(p1x-p4x);
mp1p4_ortho=-(1/mp1p4);
y_ortho1=mp1p4_ortho*(x-p1p4xCenter)+p1p4yCenter;%perndicular bisector line of point p1 and p4
y_ortho2=tan(0)*(x-a0x)+a0y;% line from A0 for getting point B0
dx=0.1;
x=-20:dx:20;
plot(x,y_ortho1);
hold on
plot(x,y_ortho2);
hold on
% Finding B0 point
eqns_function=@(var)[var(1)-(mp1p4_ortho)*(var(2)-p1p4xCenter)-p1p4yCenter;
                     var(1)-(tan(0))*(var(2)-a0x)-a0y];
initials=[0;0];
solns=fsolve(eqns_function,initials,...
             optimoptions('fsolve','Display','off'))       
b0y=solns(1);
b0x=solns(2);
plot(b0x,b0y,'r*','MarkerSize',3);
hold on
% Plotting Fixed Link A0B0
A0B0x=[a0x b0x];
A0B0y=[a0y b0y];
plot(A0B0x,A0B0y,'c','LineWidth',2);
hold on


% Locus of P1 and Locus of P4
line([p1x b0x],[p1y b0y]);
line([p4x b0x],[p4y b0y]);
shi1=atan((p1y-b0y)/(p1x-b0x));
shi2=atan((p4y-b0y)/(p4x-b0x));
shi=pi+shi1-shi2;
disp(shi);
disp(shi*180/pi);

p4locus=tan(pi-(shi/2))*(x-b0x)+b0y;
plot(x,p4locus,'--');
hold on
p1locus=tan(pi+(shi/2))*(x-b0x)+b0y;
plot(x,p1locus,'--');
hold on

% Finding A1 point
a1line=tan(-pi/3)*(x-a0x)+a0y;
plot(x,a1line,'--');
hold on
eqns_function=@(var)[var(1)-(tan(pi+(shi/2)))*(var(2)-b0x)-b0y;
                     var(1)-(tan(-pi/3))*(var(2)-a0x)-a0y];
initials=[0;0];
solns=fsolve(eqns_function,initials,...
             optimoptions('fsolve','Display','off'))       
a1y=solns(1);
a1x=solns(2);
plot(a1x,a1y,'r*','MarkerSize',3);
hold on

% Finding A4 point
a4line=tan(pi/3)*(x-a0x)+a0y;
plot(x,a4line,'--');
hold on
eqns_function=@(var)[var(1)-(tan(pi-(shi/2)))*(var(2)-b0x)-b0y;
                     var(1)-(tan(pi/3))*(var(2)-a0x)-a0y];
initials=[0;0];
solns=fsolve(eqns_function,initials,...
             optimoptions('fsolve','Display','off'))       
a4y=solns(1);
a4x=solns(2);
plot(a4x,a4y,'r*','MarkerSize',3);
hold on

% Crank A0A1 and A0A4 plotting and Crank Circle plotting
CrankLength=((a0x-a1x)^2+(a0y-a1y)^2)^0.5;
alpha=0:pi/1000:2*pi;
a0xc=CrankLength*cos(alpha)+a0x;
a0yc=CrankLength*sin(alpha)+a0y;
plot(a0xc,a0yc,'--');
hold on
line([a0x a1x],[a0y a1y],'Color','red','LineWidth',2);
hold on

% Connecting A1P1
APLength=((a1x-p1x)^2+(a1y-p1y)^2)^0.5;
A1P1x=[a1x,p1x];
A1P1y=[a1y,p1y];
plot(A1P1x,A1P1y,'r','LineWidth',2);
hold on

% Connecting A2P2 and A0A2
C1=[p2x,p2y];
C2=[a0x,a0y];
R1=APLength;
R2=CrankLength;
F=@(x) ([(x(1)-C1(1))^2+(x(2)-C1(2))^2-R1^2; ...
         (x(1)-C2(1))^2+(x(2)-C2(2))^2-R2^2]);
opt=optimoptions(@fsolve);
opt.Algorithm='levenberg-marquardt';
opt.Display='off';
x=fsolve(F,[C1(1),C1(1)+R1],opt);
fprintf('First intersection point for A2P2 line: (%f,%f)\n',x(1),x(2));
a2x=x(1);
a2y=x(2);
plot(a2x,a2y,'r*','MarkerSize',3);
x=fsolve(F,[C1(1),C1(1)-R1],opt);
fprintf('Second intersection point for A2P2 line: (%f,%f)\n',x(1),x(2));
a2x=x(1);
a2y=x(2);
plot(x(1),x(2),'r*','MarkerSize',3);
hold on

A0A2x=[a0x,a2x];
A0A2y=[a0y,a2y];
plot(A0A2x,A0A2y,'y','LineWidth',2);
hold on
A2P2x=[a2x,p2x];
A2P2y=[a2y,p2y];
plot(A2P2x,A2P2y,'y','LineWidth',2);
hold on

% Connecting A3P3 and A3A0
C1=[p3x,p3y];
C2=[a0x,a0y];
R1=APLength;
R2=CrankLength;
F=@(x) ([(x(1)-C1(1))^2+(x(2)-C1(2))^2-R1^2; ...
         (x(1)-C2(1))^2+(x(2)-C2(2))^2-R2^2]);
opt=optimoptions(@fsolve);
opt.Algorithm='levenberg-marquardt';
opt.Display='off';
x=fsolve(F,[C1(1),C1(1)+R1],opt);
fprintf('First intersection point for A3P3 line: (%f,%f)\n',x(1),x(2));
x=fsolve(F,[C1(1),C1(1)-R1],opt);
fprintf('Second intersection point for A3P3 line: (%f,%f)\n',x(1),x(2));
a3x=x(1);
a3y=x(2);
plot(x(1),x(2),'r*','MarkerSize',3);
hold on

A0A3x=[a0x,a3x];
A0A3y=[a0y,a3y];
plot(A0A3x,A0A3y,'m','LineWidth',2);
hold on
A3P3x=[a3x,p3x];
A3P3y=[a3y,p3y];
plot(A3P3x,A3P3y,'m','LineWidth',2);
hold on

% Connecting A0A4 and A4P4
A0A4x=[a0x,a4x];
A0A4y=[a0y,a4y];
plot(A0A4x,A0A4y,'g','LineWidth',2);
hold on
A4P4x=[a4x,p4x];
A4P4y=[a4y,p4y];
plot(A4P4x,A4P4y,'g','LineWidth',2);
hold on

% For finding B02 point by intersection method
P2B0Length=((p2x-b0x)^2+(p2y-b0y)^2)^0.5;   % P2A2B0 similar triangle with P1A1B02
A2B0Length=((a2x-b0x)^2+(a2y-b0y)^2)^0.5;
R1=P2B0Length;
R2=A2B0Length;
C1=[p1x,p1y];
C2=[a1x,a1y];
F=@(x) ([(x(1)-C1(1))^2+(x(2)-C1(2))^2-R1^2; ...
         (x(1)-C2(1))^2+(x(2)-C2(2))^2-R2^2]);
opt=optimoptions(@fsolve);
opt.Algorithm='levenberg-marquardt';
opt.Display='off';
x=fsolve(F,[C1(1),C1(1)+R1],opt);
fprintf('First intersection point for B02 point: (%f,%f)\n',x(1),x(2));
x=fsolve(F,[C1(1),C1(1)-R1],opt);
fprintf('Second intersection point for B02 line: (%f,%f)\n',x(1),x(2));
b02x=x(1);
b02y=x(2);
disp(b02x);
plot(b02x,b02y,'r*','MarkerSize',3);
hold on

% Finding B03 pint by intersection method
P3B0Length=((p3x-b0x)^2+(p3y-b0y)^2)^0.5;
A3B0Length=((a3x-b0x)^2+(a3y-b0y)^2)^0.5;
R1=P3B0Length;
R2=A3B0Length;
C1=[p1x,p1y];
C2=[a1x,a1y];
F=@(x) ([(x(1)-C1(1))^2+(x(2)-C1(2))^2-R1^2; ...
         (x(1)-C2(1))^2+(x(2)-C2(2))^2-R2^2]);
opt=optimoptions(@fsolve);
opt.Algorithm='levenberg-marquardt';
opt.Display='off';
x=fsolve(F,[C1(1),C1(1)+R1],opt);
fprintf('First intersection point for B03 point: (%f,%f)\n',x(1),x(2));
x=fsolve(F,[C1(1),C1(1)-R1],opt);
fprintf('Second intersection point for B03 line: (%f,%f)\n',x(1),x(2));
b03x=x(1);
b03y=x(2);
disp(b03x);
plot(b03x,b03y,'r*','MarkerSize',3);
hold on

% Connecting B0B02 and B02B03 line 
B0B02x=[b0x,b02x];
B0B02y=[b0y,b02y];
plot(B0B02x,B0B02y,'--');
hold on
B02B03x=[b02x,b03x];
B02B03y=[b02y,b03y];
plot(B02B03x,B02B03y,'--');
hold on

% Lets find the B1 point by using perpendicular bisector method
b0b02xCenter=(b0x+b02x)/2;
b0b02yCenter=(b0y+b02y)/2;
plot(b0b02xCenter,b0b02yCenter, 'r*', 'MarkerSize', 3);
hold on

b02b03xCenter=(b02x+b03x)/2;
b02b03yCenter=(b02y+b03y)/2;
plot(b02b03xCenter,b02b03yCenter, 'r*', 'MarkerSize', 3);
hold on

mb0b02=(b0y-b02y)/(b0x-b02x);
disp(mb0b02);
mb0b02_ortho=-(1/mb0b02);
disp(mb0b02_ortho);
x=-20:0.1:20;
y_ortho3=mb0b02_ortho*(x-b0b02xCenter)+b0b02yCenter;
plot(x,y_ortho3);
hold on

mb02b03=(b02y-b03y)/(b02x-b03x);
disp(mb02b03);
mb02b03_ortho=-(1/mb02b03);
disp(mb02b03_ortho);
x=-20:0.1:20;
y_ortho4=mb02b03_ortho*(x-b02b03xCenter)+b02b03yCenter;
plot(x,y_ortho4);
hold on

eqns_function=@(var)[var(1)-(mb0b02_ortho)*(var(2)-b0b02xCenter)-b0b02yCenter;
                     var(1)-(mb02b03_ortho)*(var(2)-b02b03xCenter)-b02b03yCenter];
initials=[0;0];
solns=fsolve(eqns_function,initials,...
             optimoptions('fsolve','Display','off'))       
b1y=solns(1);
b1x=solns(2);
plot(b1x,b1y,'r*','MarkerSize',3);
hold on

% Connecting Follower Link
B0B1x=[b0x,b1x];
B0B1y=[b0y,b1y];
plot(B0B1x,B0B1y,'r','LineWidth',2);
hold on

% Connecting Ternary Coupler Link
A1B1x=[a1x,b1x];
A1B1y=[a1y,b1y];
plot(A1B1x,A1B1y,'r','LineWidth',2);
hold on
P1B1x=[p1x,b1x];
P1B1y=[p1y,b1y];
plot(P1B1x,P1B1y,'r','LineWidth',2);
hold on

% Circle for follower link
FollowerLength=((b0x-b1x)^2+(b0y-b1y)^2)^0.5;
alpha=0:pi/1000:2*pi;
b0xc=FollowerLength*cos(alpha)+b0x;
b0yc=FollowerLength*sin(alpha)+b0y;
plot(b0xc,b0yc,'--');
hold on

% Connecting P2B2 and A2B2
P1B1Length=((p1x-b1x)^2+(p1y-b1y)^2)^0.5;
C1=[p2x,p2y];
C2=[b0x,b0y];
R1=P1B1Length;
R2=FollowerLength;
F=@(x) ([(x(1)-C1(1))^2+(x(2)-C1(2))^2-R1^2; ...
         (x(1)-C2(1))^2+(x(2)-C2(2))^2-R2^2]);
opt=optimoptions(@fsolve);
opt.Algorithm='levenberg-marquardt';
opt.Display='off';
x=fsolve(F,[C1(1),C1(1)+R1],opt);
fprintf('First intersection point for P2B2 line: (%f,%f)\n',x(1),x(2));
b2x=x(1);
b2y=x(2);
plot(b2x,b2y,'r*','MarkerSize',3);
hold on
x=fsolve(F,[C1(1),C1(1)-R1],opt);
fprintf('Second intersection point for P2B2 line: (%f,%f)\n',x(1),x(2));

B0B2x=[b0x,b2x];
B0B2y=[b0y,b2y];
plot(B0B2x,B0B2y,'y','LineWidth',2);
hold on
A2B2x=[a2x,b2x];
A2B2y=[a2y,b2y];
plot(A2B2x,A2B2y,'y','LineWidth',2);
hold on
P2B2x=[p2x,b2x];
P2B2y=[p2y,b2y];
plot(P2B2x,P2B2y,'y','LineWidth',2);
hold on

% Connecting P3B3 and A3B3
C1=[p3x,p3y];
C2=[b0x,b0y];
R1=P1B1Length;
R2=FollowerLength;
F=@(x) ([(x(1)-C1(1))^2+(x(2)-C1(2))^2-R1^2; ...
         (x(1)-C2(1))^2+(x(2)-C2(2))^2-R2^2]);
opt=optimoptions(@fsolve);
opt.Algorithm='levenberg-marquardt';
opt.Display='off';
x=fsolve(F,[C1(1),C1(1)+R1],opt);
fprintf('First intersection point for P2B2 line: (%f,%f)\n',x(1),x(2));
b3x=x(1);
b3y=x(2);
plot(b3x,b3y,'r*','MarkerSize',3);
hold on
x=fsolve(F,[C1(1),C1(1)-R1],opt);
fprintf('Second intersection point for P2B2 line: (%f,%f)\n',x(1),x(2));

B0B3x=[b0x,b3x];
B0B3y=[b0y,b3y];
plot(B0B3x,B0B3y,'m','LineWidth',2);
hold on
A3B3x=[a3x,b3x];
A3B3y=[a3y,b3y];
plot(A3B3x,A3B3y,'m','LineWidth',2);
hold on
P3B3x=[p3x,b3x];
P3B3y=[p3y,b3y];
plot(P3B3x,P3B3y,'m','LineWidth',2);
hold on

% Connecting P4B4 and A4B4
C1=[p4x,p4y];
C2=[b0x,b0y];
R1=P1B1Length;
R2=FollowerLength;
F=@(x) ([(x(1)-C1(1))^2+(x(2)-C1(2))^2-R1^2; ...
         (x(1)-C2(1))^2+(x(2)-C2(2))^2-R2^2]);
opt=optimoptions(@fsolve);
opt.Algorithm='levenberg-marquardt';
opt.Display='off';
x=fsolve(F,[C1(1),C1(1)+R1],opt);
fprintf('First intersection point for P2B2 line: (%f,%f)\n',x(1),x(2));
b4x=x(1);
b4y=x(2);
plot(b4x,b4y,'r*','MarkerSize',3);
hold on
x=fsolve(F,[C1(1),C1(1)-R1],opt);
fprintf('Second intersection point for P2B2 line: (%f,%f)\n',x(1),x(2));

B0B4x=[b0x,b4x];
B0B4y=[b0y,b4y];
plot(B0B4x,B0B4y,'g','LineWidth',2);
hold on
A4B4x=[a4x,b4x];
A4B4y=[a4y,b4y];
plot(A4B4x,A4B4y,'g','LineWidth',2);
hold on
P4B4x=[p4x,b4x];
P4B4y=[p4y,b4y];
plot(P4B4x,P4B4y,'g','LineWidth',2);
hold on

A0A1Length=((a0x-a1x)^2+(a0y-a1y)^2)^0.5;
A1P1Length=((a1x-p1x)^2+(a1y-p1y)^2)^0.5;
P1B1Length=((p1x-b1x)^2+(p1y-b1y)^2)^0.5;
A1B1Length=((a1x-b1x)^2+(a1y-b1y)^2)^0.5;
B1B0Length=((b1x-b0x)^2+(b1y-b0y)^2)^0.5;
fprintf('A0A1=%f\t\tA1P1=%f\t\tP1B1=%f\t\tA1B1=%f\t\tB1B0=%f\n',A0A1Length,A1P1Length,P1B1Length,A1B1Length,B1B0Length);

A0A2Length=((a0x-a2x)^2+(a0y-a2y)^2)^0.5;
A2P2Length=((a2x-p2x)^2+(a2y-p2y)^2)^0.5;
P2B2Length=((p2x-b2x)^2+(p2y-b2y)^2)^0.5;
A2B2Length=((a2x-b2x)^2+(a2y-b2y)^2)^0.5;
B2B0Length=((b2x-b0x)^2+(b2y-b0y)^2)^0.5;
fprintf('A0A2=%f\t\tA2P2=%f\t\tP2B2=%f\t\tA2B2=%f\t\tB2B0=%f\n',A0A2Length,A2P2Length,P2B2Length,A2B2Length,B2B0Length);

A0A3Length=((a0x-a3x)^2+(a0y-a3y)^2)^0.5;
A3P3Length=((a3x-p3x)^2+(a3y-p3y)^2)^0.5;
P3B3Length=((p3x-b3x)^2+(p3y-b3y)^2)^0.5;
A3B3Length=((a3x-b3x)^2+(a3y-b3y)^2)^0.5;
B3B0Length=((b3x-b0x)^2+(b3y-b0y)^2)^0.5;
fprintf('A0A3=%f\t\tA3P3=%f\t\tP3B3=%f\t\tA3B3=%f\t\tB3B0=%f\n',A0A3Length,A3P3Length,P3B3Length,A3B3Length,B3B0Length);

A0A4Length=((a0x-a4x)^2+(a0y-a4y)^2)^0.5;
A4P4Length=((a4x-p4x)^2+(a4y-p4y)^2)^0.5;
P4B4Length=((p4x-b4x)^2+(p4y-b4y)^2)^0.5;
A4B4Length=((a4x-b4x)^2+(a4y-b4y)^2)^0.5;
B4B0Length=((b4x-b0x)^2+(b4y-b0y)^2)^0.5;
fprintf('A0A4=%f\t\tA4P4=%f\t\tP4B4=%f\t\tA4B4=%f\t\tB4B0=%f\n',A0A4Length,A4P4Length,P4B4Length,A4B4Length,B4B0Length);
disp(b0x);
disp(b0y);
disp(a1x);
disp(a1y);
disp(b1x);
disp(b1y);

d=b0x-a0x;
a=CrankLength;
b=A1B1Length;
c=FollowerLength;
A0=[a0x a0y];
B0=[b0x b0y];
k=0.5;
for t=1:500
    theta2=k*(t/10);
    A1=A0+a*[cos(theta2) sin(theta2)];
    k1=d/a;
    k2=d/c;
    k3=(a*a-b*b+c*c+d*d)/(2*a*c);
    k4=d/b;
    k5=(c*c-d*d-a*a-b*b)/(2*a*b);
    A=cos(theta2)-k1-k2*cos(theta2)+k3;
    B=-2*sin(theta2);
    C=k1-(k2+1)*cos(theta2)+k3;
    D=cos(theta2)-k1+k4*cos(theta2)+k5;
    E=-2*sin(theta2);
    F=k1+(k4-1)*cos(theta2)+k5;
    theta3=2*atan(real((-E-(E*E-4*D*F)^0.5)/2*D));
    theta4=2*atan(real((-B-(B*B-4*A*C)^0.5)/2*A));
    B1=B0+c*[cos(theta4) sin(theta4)];  
    
    % Getting P1
    C1=[A1(1),A1(2)];
    C2=[B1(1),B1(2)];
    R1=A1P1Length;
    R2=P1B1Length;
    F=@(x) ([(x(1)-C1(1))^2+(x(2)-C1(2))^2-R1^2; ...
           (x(1)-C2(1))^2+(x(2)-C2(2))^2-R2^2]);
    opt=optimoptions(@fsolve);
    opt.Algorithm='levenberg-marquardt';
    opt.Display='off';
    x=fsolve(F,[C1(1),C1(1)+R1],opt);
    fprintf('First intersection point for P2B2 line: (%f,%f)\n',x(1),x(2));
    P1(1)=x(1);
    P1(2)=x(2);
    plot(b3x,b3y,'r*','MarkerSize',3);
    hold on
    x=fsolve(F,[C1(1),C1(1)-R1],opt);
    fprintf('Second intersection point for P2B2 line: (%f,%f)\n',x(1),x(2));
    P1=[P1(1) P1(2)];
    
    crank=line([A0(1) A1(1)],[A0(2) A1(2)]);
    coupler=line([A1(1) B1(1)],[A1(2) B1(2)]);
    follower=line([B0(1) B1(1)],[B0(2) B1(2)]);
    coupler2=line([A1(1) P1(1)],[A1(2) P1(2)]);
    coupler3=line([B1(1) P1(1)],[B1(2) P1(2)]);
    
    A1_traj=viscircles([a0x a0y],a,'LineStyle','--');
    A0_circ=viscircles(A0,0.1);
    A1_circ=viscircles(A1,0.1);
    B1_traj=viscircles([b0x b0y],c,'LineStyle','--');
    B1_circ=viscircles(B1,0.1);
    B0_circ=viscircles(B0,0.1);
    P1_circ=viscircles(P1,0.1);
    
    pause(0.001);
    delete(crank);
    delete(coupler);
    delete(follower);
    delete(coupler2);
    delete(coupler3);
    delete(A0_circ);
    delete(A1_circ);
    delete(B1_circ);
    delete(B0_circ);
    delete(P1_circ);
    
    %Check length
    cranklength=((A0(1)-A1(1))^2+(A0(2)-A1(2))^2)^0.5;
    couplerlength=((B1(1)-A1(1))^2+(B1(2)-A1(2))^2)^0.5;
    followerlength=((B0(1)-B1(1))^2+(B0(2)-B1(2))^2)^0.5;
    coupler2length=((P1(1)-A1(1))^2+(P1(2)-A1(2))^2)^0.5;
    coupler3length=((B1(1)-P1(1))^2+(B1(2)-P1(2))^2)^0.5;
    
    disp(cranklength);
    disp(couplerlength);
    disp(followerlength);
    disp(coupler2length);
    disp(coupler3length);
end




axis([-10 10 -10 10])
axis('equal')


