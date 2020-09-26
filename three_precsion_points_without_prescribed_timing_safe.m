 %%%%%%%%% Path Generation Three Precision Point Problem without Prescribed Timing %%%%%%%%%%
% Fixed line A0B0
a0x=-2;%a0x=input('Enter x coordinate of A0 point=');
a0y=2;%a0y=input('Enter y coordinate of A0 point=');
b0x=7;%b0x=input('Enter x coordinate of B0 point=');
b0y=2;%b0y=input('Enter y coordinate of B0 point=');

% Different Position
p1x=3;%p1x=input('Enter x coordinate of P1 point=');
p1y=8;%p1y=input('Enter y coordinate of P1 point=');

p2x=4;%p2x=input('Enter x coordinate of P2 point=');
p2y=8.5;%p2y=input('Enter y coordinate of P2 point=');

p3x=5;%p3x=input('Enter x coordinate of P3 point=');
p3y=8.2;%p3y=input('Enter y coordinate of P3 point=');

a1x=-4;%p3x=input('Enter x coordinate of A1 point=');
a1y=3;%p3y=input('Enter y coordinate of A1 point=');

% Connect different links
A0B0x=[a0x,b0x];
A0B0y=[a0y,b0y];

% Mechanism Ploting
plot(A0B0x,A0B0y,'green','LineWidth',2);
hold on

% Crank and Crank circle Plotting
A0A1x=[a0x,a1x];
A0A1y=[a0y,a1y];
CrankLength=((a0x-a1x)^2+(a0y-a1y)^2)^0.5;
alpha=0:pi/1000:2*pi;
a0xc=CrankLength*cos(alpha)+a0x;
a0yc=CrankLength*sin(alpha)+a0y;
plot(a0xc,a0yc,'--');
hold on

% Connecting A1P1
APLength=((a1x-p1x)^2+(a1y-p1y)^2)^0.5;
plot(A0A1x,A0A1y,'r','LineWidth',2);
hold on
A1P1x=[a1x,p1x];
A1P1y=[a1y,p1y];
plot(A1P1x,A1P1y,'r','LineWidth',2);
hold on

% Connecting A2P2
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

% Connecting A3P3
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
a3x=x(1);
a3y=x(2);
plot(a3x,a3y,'r*','MarkerSize',3);
hold on
x=fsolve(F,[C1(1),C1(1)-R1],opt);
fprintf('Second intersection point for A3P3 line: (%f,%f)\n',x(1),x(2));
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
x=-8:0.1:15;
y_ortho3=mb0b02_ortho*(x-b0b02xCenter)+b0b02yCenter;
plot(x,y_ortho3);
hold on

mb02b03=(b02y-b03y)/(b02x-b03x);
disp(mb02b03);
mb02b03_ortho=-(1/mb02b03);
disp(mb02b03_ortho);
x=-8:0.1:15;
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
A1B1x=[a1x,b1x];
A1B1y=[a1y,b1y];
plot(A1B1x,A1B1y,'r','LineWidth',2);
hold on
P1B1x=[p1x,b1x];
P1B1y=[p1y,b1y];
plot(P1B1x,P1B1y,'r','LineWidth',2);
hold on

FollowerLength=((b0x-b1x)^2+(b0y-b1y)^2)^0.5;
% Circle for follower link
alpha=0:pi/1000:2*pi;
b0xc=FollowerLength*cos(alpha)+b0x;
b0yc=FollowerLength*sin(alpha)+b0y;
plot(b0xc,b0yc,'--');
hold on

% Length of P1B1
P1B1Length=((p1x-b1x)^2+(p1y-b1y)^2)^0.5;

% Connecting P2B2
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
plot(x(1),x(2),'r*','MarkerSize',3);
hold on

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

% Connecting P3B3
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
plot(x(1),x(2),'r*','MarkerSize',3);
hold on

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

A1P1Length=((a1x-p1x)^2+(a1y-p1y)^2)^0.5;
P1B1Length=((p1x-b1x)^2+(p1y-b1y)^2)^0.5;
A1B1Length=((a1x-b1x)^2+(a1y-b1y)^2)^0.5;
fprintf('A1P1Length=%f\t\tP1B1Length=%f\t\tA1B1Length=%f\n',A1P1Length,P1B1Length,A1B1Length);

A2P2Length=((a2x-p2x)^2+(a2y-p2y)^2)^0.5;
P2B2Length=((p2x-b2x)^2+(p2y-b2y)^2)^0.5;
A2B2Length=((a2x-b2x)^2+(a2y-b2y)^2)^0.5;
fprintf('A2P2Length=%f\t\tP2B2Length=%f\t\tA2B2Length=%f\n',A2P2Length,P2B2Length,A2B2Length);

A3P3Length=((a3x-p3x)^2+(a3y-p3y)^2)^0.5;
P3B3Length=((p3x-b3x)^2+(p3y-b3y)^2)^0.5;
A3B3Length=((a3x-b3x)^2+(a3y-b3y)^2)^0.5;
fprintf('A3P3Length=%f\t\tP3B3Length=%f\t\tA3B3Length=%f\n',A3P3Length,P3B3Length,A3B3Length);


d=b0x-a0x;
a=CrankLength;
b=A1B1Length;
c=FollowerLength;
A0=[a0x a0y];
B0=[b0x b0y];

k=0.5;
for t=1:360
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
    theta3=2*atan(real((-E+(E*E-4*D*F)^0.5)/2*D));
    theta4=2*atan(real((-B+(B*B-4*A*C)^0.5)/2*A));
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
end

axis([-20 20 -20 20])
axis('equal')
