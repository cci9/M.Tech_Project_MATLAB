%%%%%%%%% Path Generation Three Precision Point Problem with Prescribed Timing %%%%%%%%%%
% Fixed line A0B0
% Different Position
p1x=input('Enter x coordinate of P1 point=');
p1y=input('Enter y coordinate of P1 point=');

p2x=input('Enter x coordinate of P2 point=');
p2y=input('Enter y coordinate of P2 point=');

p3x=input('Enter x coordinate of P3 point=');
p3y=input('Enter y coordinate of P3 point=');

p4x=input('Enter x coordinate of P4 point=');
p4y=input('Enter y coordinate of P4 point=');

% Rotational angle of crank with respect to different position
theta12=input('Enter crank rotation angle in degrees when coupler moves from point P1 to point P2=');
thetal21=thetal12;
theta21=thetal21*pi/180; %%degree to radian angle 

theta13=input('Enter crank rotation angle in degrees when coupler moves from point P1 to point P3=');
thetal31=thetal13;
theta31=thetal31*pi/180; %% degree to radian angle 

theta14=input('Enter crank rotation angle in degrees when coupler moves from point P1 to point P4=');
thetal41=thetal14;
theta41=thetal41*pi/180; %% degree to radian angle

% Lets plot the point through which coupler sgould pass
plot(p1x,p1y, 'r*', 'MarkerSize', 3);
hold on
plot(p2x,p2y, 'r*', 'MarkerSize', 3);
hold on
plot(p3x,p3y, 'r*', 'MarkerSize', 3);
hold on
plot(p4x,p4y, 'r*', 'MarkerSize', 3);
hold on

% Connect position of P2 with P1
P1xP2x=[p1x,p2x];
P1yP2y=[p1y,p2y];
plot(P1xP2x,P1yP2y,'--');
hold on

% Finding perpendicular bisector of P1P2
p1p2xCenter=(p1x+p2x)/2;
p1p2yCenter=(p1y+p2y)/2;
plot(p1p2xCenter, p1p2yCenter, 'r*', 'MarkerSize', 3);
hold on

mp1p2=(p1y-p2y)/(p1x-p2x);
mp1p2_ortho=-(1/mp1p2);
y_ortho_p1p2=mp1p2_ortho*(x-p1p2xCenter)+p1p2yCenter;

% Lets find A0 point on C1C2 perpendicular bisector line
y_p1p2=mp1p2*(x-p1x)+p1y;
theta_p1p2=atan(mp1p2);
theta_p1a0=theta_p1p2-(pi-theta21)/2;
y_p1a0=(tan(theta_p1a0))*(x-p1x)+p1y;
disp(theta_p1p2);
disp(theta_p1a0);

dx=0.1;
x=-300:dx:300;
plot(x,y_ortho_p1p2,'--');
hold on
plot(x,y_p1p2,'--');
hold on
plot(x,y_p1a0,'--');
hold on

% Lets find A0 point by intersection of  P1P2 perp. bisec. and P1A0
y_ortho_p1p2=mp1p2_ortho*(x-p1p2xCenter)+p1p2yCenter;
y_p1a0=(tan(theta_p1a0))*(x-p1x)+p1y;
eqns_function=@(var)[var(1)-(mp1p2_ortho)*(var(2)-p1p2xCenter)-p1p2yCenter;
                     var(1)-(tan(theta_p1a0))*(var(2)-p1x)-p1y];
initials=[0;0];
solns=fsolve(eqns_function,initials,...
             optimoptions('fsolve','Display','off'))       
a0y=solns(1);
a0x=solns(2);
plot(a0x,a0y,'r*','MarkerSize',15);
hold on

% Lets plot A0P1, A0P2, A0P3, and A0P4
A0P1x=[a0x,p1x];
A0P1y=[a0y,p1y];
A0P2x=[a0x,p2x];
A0P2y=[a0y,p2y];
A0P3x=[a0x,p3x];
A0P3y=[a0y,p3y];
A0P4x=[a0x,p4x];
A0P4y=[a0y,p4y];
plot(A0P1x,A0P1y,'--');
hold on
plot(A0P2x,A0P2y,'--');
hold on
plot(A0P3x,A0P3y,'--');
hold on
plot(A0P4x,A0P4y,'--');
hold on

% Find the length and angle of each position point with fixed point
L3=((p3x-a0x)^2+(p3y-a0y)^2)^0.5;
disp(L3);
alphaa0p3=atan((p3y-a0y)/(p3x-a0x));

L4=((p4x-a0x)^2+(p4y-a0y)^2)^0.5;
disp(L4);
alphaa0p4=atan((p4y-a0y)/(p4x-a0x));

% Find the different poition of P2, P3 and P4 after rotation by thetal21,thetal31, and thetal41
p2nx=p1x;
p2ny=p1y;

p3nx=a0x+L3*cos(theta31+alphaa0p3); 
p3ny=a0y+L3*sin(theta31+alphaa0p3);
A0P3nx=[a0x,p3nx];
A0P3ny=[a0y,p3ny];
plot(p3nx,p3ny,'r*','Markersize',3);
hold on
plot(A0P3nx,A0P3ny,'--');
hold on
L3n=((p3nx-a0x)^2+(p3ny-a0y)^2)^0.5;
disp(L3n);

p4nx=a0x+L4*cos(theta41+alphaa0p4);
p4ny=a0y+L4*sin(theta41+alphaa0p4);
A0P4nx=[a0x,p4nx];
A0P4ny=[a0y,p4ny];
plot(p4nx,p4ny,'r*','Markersize',3);
hold on
plot(A0P4nx,A0P4ny,'--');
hold on
L4n=((p4nx-a0x)^2+(p4ny-a0y)^2)^0.5;
disp(L4n);

% Connect new posi. of P3 with P1 and new posi. of P4 with new posi. of P3
P1xP3nx=[p1x,p3nx];
P1yP3ny=[p1y,p3ny];
plot(P1xP3nx,P1yP3ny,'--');
hold on
P3nxP4nx=[p3nx,p4nx];
P3nyP4ny=[p3ny,p4ny];
plot(P3nxP4nx,P3nyP4ny,'--');
hold on

% Lets find perpendicular bisector of P1P3n and P3nP4n
p1p3nxCenter=(p1x+p3nx)/2;
p1p3nyCenter=(p1y+p3ny)/2;
plot(p1p3nxCenter, p1p3nyCenter, 'r*', 'MarkerSize', 3);
hold on

p3np4nxCenter=(p3nx+p4nx)/2;
p3np4nyCenter=(p3ny+p4ny)/2;
plot(p3np4nxCenter, p3np4nyCenter, 'r*', 'MarkerSize', 3);
hold on

mp1p3n=(p1y-p3ny)/(p1x-p3nx);
mp1p3n_ortho=-(1/mp1p3n);
y_ortho_p1p3n=mp1p3n_ortho*(x-p1p3nxCenter)+p1p3nyCenter;

mp3np4n=(p3ny-p4ny)/(p3nx-p4nx);
mp3np4n_ortho=-(1/mp3np4n);
y_ortho_p3np4n=mp3np4n_ortho*(x-p3np4nxCenter)+p3np4nyCenter;
dx=0.1;
x=-300:dx:300;
plot(x,y_ortho_p1p3n,'--');
hold on
plot(x,y_ortho_p3np4n,'--');
hold on

% Lets find intersection point between P1P3n and P3nP4n line
y_ortho_p1p3n=mp1p3n_ortho*(x-p1p3nxCenter)+p1p3nyCenter;
y_ortho_p3np4n=mp3np4n_ortho*(x-p3np4nxCenter)+p3np4nyCenter;
eqns_function=@(var)[var(1)-(mp1p3n_ortho)*(var(2)-p1p3nxCenter)-p1p3nyCenter;
                     var(1)-(mp3np4n_ortho)*(var(2)-p3np4nxCenter)-p3np4nyCenter];
initials=[0;0];
solns=fsolve(eqns_function,initials,...
             optimoptions('fsolve','Display','off'))       
a1y=solns(1);
a1x=solns(2);
plot(a1x,a1y,'r*','MarkerSize',3);
hold on

% Crank and Crank circle Plotting
CrankLength=((a0x-a1x)^2+(a0y-a1y)^2)^0.5;
A0A1x=[a0x,a1x];
A0A1y=[a0y,a1y];
alpha=0:pi/1000:2*pi;
a0xc=CrankLength*cos(alpha)+a0x;
a0yc=CrankLength*sin(alpha)+a0y;
plot(a0xc,a0yc,'--');
hold on

% Connecting A0A1 and A1P1
plot(A0A1x,A0A1y,'r','LineWidth',2);
hold on
A1P1x=[a1x,p1x];
A1P1y=[a1y,p1y];
plot(A1P1x,A1P1y,'r','LineWidth',2);
hold on


% Connecting A2P2
A1P1Length=((a1x-p1x)^2+(a1y-p1y)^2)^0.5;
A0A1Length=((a0x-a1x)^2+(a0y-a1y)^2)^0.5;
% Center of the circles
C1=[a0x;a0y];
C2=[p2x;p2y];
% Radius of the circles
R1=A0A1Length;
R2=A1P1Length;
d2 = sum((C2-C1).^2);
   C0 = (C1+C2)/2+(R1^2-R2^2)/d2/2*(C2-C1);
   t = ((R1+R2)^2-d2)*(d2-(R2-R1)^2);
   if t <= 0
     fprintf('The circles don''t intersect.\n')
   else
     T = (sqrt(t)/d2/2)*[0 -1;1 0]*(C2-C1);
     Ca = C0 + T; % Pa and Pb are circles' intersection points
     Cb = C0 - T;
   end
   disp(Ca);
   disp(Cb);
a2x1=Ca(1);
a2y1=Ca(2);
a2x2=Cb(1);
a2y2=Cb(2);
plot(a2x1,a2y1,'r*','MarkerSize',3);
hold on
plot(a2x2,a2y2,'r*','MarkerSize',3);
hold on
% Checking correct A2 point
A1A2distance=((a1x-a2x2)^2+(a1y-a2y2)^2)^0.5;
alph21=acos((2*CrankLength^2-A1A2distance^2)/(2*CrankLength^2));
disp(alph21);
disp(theta21);
if alph21==theta21
    a2x=a2x2;
    a2y=a2y2;
else
    a2x=a2x1;
    a2y=a2y1;
end
disp([a2x a2y]);
% Connecting A0A2 and A2P2
A0A2x=[a0x,a2x];
A0A2y=[a0y,a2y];
plot(A0A2x,A0A2y,'y','LineWidth',2);
hold on
A2P2x=[a2x,p2x];
A2P2y=[a2y,p2y];
plot(A2P2x,A2P2y,'y','LineWidth',2);
hold on

% Connecting A3P3
A1P1Length=((a1x-p1x)^2+(a1y-p1y)^2)^0.5;
A0A1Length=((a0x-a1x)^2+(a0y-a1y)^2)^0.5;
% Center of the circles
C1=[a0x;a0y];
C2=[p3x;p3y];
% Radius of the circles
R1=A0A1Length;
R2=A1P1Length;
d2 = sum((C2-C1).^2);
   C0 = (C1+C2)/2+(R1^2-R2^2)/d2/2*(C2-C1);
   t = ((R1+R2)^2-d2)*(d2-(R2-R1)^2);
   if t <= 0
     fprintf('The circles don''t intersect.\n')
   else
     T = (sqrt(t)/d2/2)*[0 -1;1 0]*(C2-C1);
     Ca = C0 + T; % Ca and Cb are circles' intersection points
     Cb = C0 - T;
   end
   disp(Ca);
   disp(Cb);
a3x1=Ca(1);
a3y1=Ca(2);
a3x2=Cb(1);
a3y2=Cb(2);
disp([a3x1 a3y1]);
disp([a3x2 a3y2]);
plot(a3x1,a3y1,'r*','MarkerSize',3);
hold on
plot(a3x2,a3y2,'r*','MarkerSize',3);
hold on
% Checking correct A3 point
A1A3distance=((a1x-a3x2)^2+(a1y-a3y2)^2)^0.5;
alph31=acos((2*CrankLength^2-A1A2distance^2)/(2*CrankLength^2));
disp(alph31);
disp(theta31);
if alph31==theta31
    a3x=a3x2;
    a3y=a3y2;
else
    a3x=a3x1;
    a3y=a3y1;
end
disp([a3x a3y]);
% Connecting A0A3 and A3P3
A0A3x=[a0x,a3x];
A0A3y=[a0y,a3y];
plot(A0A3x,A0A3y,'m','LineWidth',2);
hold on
A3P3x=[a3x,p3x];
A3P3y=[a3y,p3y];
plot(A3P3x,A3P3y,'m','LineWidth',2);
hold on

% Connecting A4P4
A1P1Length=((a1x-p1x)^2+(a1y-p1y)^2)^0.5;
A0A1Length=((a0x-a1x)^2+(a0y-a1y)^2)^0.5;
% Center of the circles
C1=[a0x;a0y];
C2=[p4x;p4y];
% Radius of the circles
R1=A0A1Length;
R2=A1P1Length;
d2 = sum((C2-C1).^2);
   C0 = (C1+C2)/2+(R1^2-R2^2)/d2/2*(C2-C1);
   t = ((R1+R2)^2-d2)*(d2-(R2-R1)^2);
   if t <= 0
     fprintf('The circles don''t intersect.\n')
   else
     T = (sqrt(t)/d2/2)*[0 -1;1 0]*(C2-C1);
     Ca = C0 + T; % Ca and Cb are circles' intersection points
     Cb = C0 - T;
   end
   disp(Ca);
   disp(Cb);
a4x1=Ca(1);
a4y1=Ca(2);
a4x2=Cb(1);
a4y2=Cb(2);
disp([a4x1 a4y1]);
disp([a4x2 a4y2]);
plot(a4x1,a4y1,'r*','MarkerSize',3);
hold on
plot(a4x2,a4y2,'r*','MarkerSize',3);
hold on
% Checking correct A4 point
A1A4distance=((a1x-a4x2)^2+(a1y-a4y2)^2)^0.5;
alph41=acos((2*CrankLength^2-A1A2distance^2)/(2*CrankLength^2));
disp(alph41);
disp(theta41);
if alph41==theta41
    a4x=a4x2;
    a4y=a4y2;
else
    a4x=a4x1;
    a4y=a4y1;
end
disp([a4x a4y]);
% Connecting A0A4 and A4P4
A0A4x=[a0x,a4x];
A0A4y=[a0y,a4y];
plot(A0A4x,A0A4y,'g','LineWidth',2);
hold on
A4P4x=[a4x,p4x];
A4P4y=[a4y,p4y];
plot(A4P4x,A4P4y,'g','LineWidth',2);
hold on

% Lets find perpendicular bisector of P3P4 and A3A4
P3P4x=[p3x,p4x];
P3P4y=[p3y,p4y];
plot(P3P4x,P3P4y,'--');
hold on
p3p4xCenter=(p3x+p4x)/2;
p3p4yCenter=(p3y+p4y)/2;
plot(p3p4xCenter, p3p4yCenter, 'r*', 'MarkerSize', 3);
hold on

A3A4x=[a3x,a4x];
A3A4y=[a3y,a4y];
plot(A3A4x,A3A4y,'--');
hold on
a3a4xCenter=(a3x+a4x)/2;
a3a4yCenter=(a3y+a4y)/2;
plot(a3a4xCenter, a3a4yCenter, 'r*', 'MarkerSize', 3);
hold on

mp3p4=(p4y-p3y)/(p4x-p3x);
mp3p4_ortho=-(1/mp3p4);
y_ortho_p3p4=(mp3p4_ortho)*(x-p3p4xCenter)+p3p4yCenter;
plot(x,y_ortho_p3p4,'--');
hold on

ma3a4=(a4y-a3y)/(a4x-a3x);
ma3a4_ortho=-(1/ma3a4);
y_ortho_a3a4=(ma3a4_ortho)*(x-a3a4xCenter)+a3a4yCenter;
plot(x,y_ortho_a3a4,'--');
hold on


% Lets find intersection point between P3P4 and A3A4 line
y_ortho_p3p4=(mp3p4_ortho)*(x-p3p4xCenter)+p3p4yCenter;
y_ortho_a3a4=(ma3a4_ortho)*(x-a3a4xCenter)+a3a4yCenter;
eqns_function=@(var)[var(1)-(mp3p4_ortho)*(var(2)-p3p4xCenter)-p3p4yCenter;
                     var(1)-(ma3a4_ortho)*(var(2)-a3a4xCenter)-a3a4yCenter];
initials=[0;0];
solns=fsolve(eqns_function,initials,...
             optimoptions('fsolve','Display','off'))       
b3y=solns(1);
b3x=solns(2);
plot(b3x,b3y,'r*','MarkerSize',3);
hold on
disp([b3x b3y]);

% Lets plot triangle A3B3P3
A3B3x=[a3x,b3x];
A3B3y=[a3y,b3y];
plot(A3B3x,A3B3y,'m','LineWidth',2);
hold on
P3B3x=[p3x,b3x];
P3B3y=[p3y,b3y];
plot(P3B3x,P3B3y,'m','LineWidth',2);
hold on

% Lets plot triangle A4B4P4
b4x=b3x;
b4y=b3y;
A4B4x=[a4x,b4x];
A4B4y=[a4y,b4y];
plot(A4B4x,A4B4y,'g','LineWidth',2);
hold on
P4B4x=[p4x,b4x];
P4B4y=[p4y,b4y];
plot(P4B4x,P4B4y,'g','LineWidth',2);
hold on

% Lets find the point B2 and plot triangle A2B2C2
A4B4Length=((a4x-b4x)^2+(a4y-b4y)^2)^0.5;
P4B4Length=((p4x-b4x)^2+(p4y-b4y)^2)^0.5;
% Center of the circles
C1=[a2x;a2y];
C2=[p2x;p2y];
% Radius of the circles
R1=A4B4Length;
R2=P4B4Length;
d2 = sum((C2-C1).^2);
   C0 = (C1+C2)/2+(R1^2-R2^2)/d2/2*(C2-C1);
   t = ((R1+R2)^2-d2)*(d2-(R2-R1)^2);
   if t <= 0
     fprintf('The circles don''t intersect.\n')
   else
     T = (sqrt(t)/d2/2)*[0 -1;1 0]*(C2-C1);
     Ca = C0 + T; % Ca and Cb are circles' intersection points
     Cb = C0 - T;
   end
   disp(Ca);
   disp(Cb);
b2x=Ca(1);
b2y=Ca(2);
disp([b2x b2y]);

A2B2x=[a2x,b2x];
A2B2y=[a2y,b2y];
plot(A2B2x,A2B2y,'y','LineWidth',2);
hold on
P2B2x=[p2x,b2x];
P2B2y=[p2y,b2y];
plot(P2B2x,P2B2y,'y','LineWidth',2);
hold on

% Lets find the point B1 and plot triangle A1B1C1
A4B4Length=((a4x-b4x)^2+(a4y-b4y)^2)^0.5;
P4B4Length=((p4x-b4x)^2+(p4y-b4y)^2)^0.5;
% Center of the circles
C1=[a1x;a1y];
C2=[p1x;p1y];
% Radius of the circles
R1=A4B4Length;
R2=P4B4Length;
d2 = sum((C2-C1).^2);
   C0 = (C1+C2)/2+(R1^2-R2^2)/d2/2*(C2-C1);
   t = ((R1+R2)^2-d2)*(d2-(R2-R1)^2);
   if t <= 0
     fprintf('The circles don''t intersect.\n')
   else
     T = (sqrt(t)/d2/2)*[0 -1;1 0]*(C2-C1);
     Ca = C0 + T; % Ca and Cb are circles' intersection points
     Cb = C0 - T;
   end
   disp(Ca);
   disp(Cb);
b1x=Ca(1);
b1y=Ca(2);
disp([b1x b1y]);

A1B1x=[a1x,b1x];
A1B1y=[a1y,b1y];
plot(A1B1x,A1B1y,'r','LineWidth',2);
hold on
P1B1x=[p1x,b1x];
P1B1y=[p1y,b1y];
plot(P1B1x,P1B1y,'r','LineWidth',2);
hold on

% Lets find perpendicular bisector of B1B2 and B2B3
B1B2x=[b1x,b2x];
B1B2y=[b1y,b2y];
plot(B1B2x,B1B2y,'--');
hold on
b1b2xCenter=(b1x+b2x)/2;
b1b2yCenter=(b1y+b2y)/2;
plot(b1b2xCenter,b1b2yCenter, 'r*', 'Markersize',3);
hold on

B3B2x=[b3x,b2x];
B3B2y=[b3y,b2y];
plot(B3B2x,B3B2y,'--');
hold on
b3b2xCenter=(b3x+b2x)/2;
b3b2yCenter=(b3y+b2y)/2;
plot(b3b2xCenter,b3b2yCenter, 'r*', 'Markersize',3);
hold on

mb1b2=(b1y-b2y)/(b1x-b2x);
mb1b2_ortho=-(1/mb1b2);
y_ortho_b1b2=mb1b2_ortho*(x-b1b2xCenter)+b1b2yCenter;
plot(x,y_ortho_b1b2,'--');
hold on

mb3b2=(b3y-b2y)/(b3x-b2x);
mb3b2_ortho=-(1/mb3b2);
y_ortho_b3b2=mb3b2_ortho*(x-b3b2xCenter)+b3b2yCenter;
plot(x,y_ortho_b3b2,'--');
hold on

y_ortho_b1b2=mb1b2_ortho*(x-b1b2xCenter)+b1b2yCenter;
y_ortho_b3b2=mb3b2_ortho*(x-b3b2xCenter)+b3b2yCenter;
eqns_function=@(var)[var(1)-(mb1b2_ortho)*(var(2)-b1b2xCenter)-b1b2yCenter;
                     var(1)-(mb3b2_ortho)*(var(2)-b3b2xCenter)-b3b2yCenter];
initials=[0;0];
solns=fsolve(eqns_function,initials,...
             optimoptions('fsolve','Display','off'))       
b0y=solns(1);
b0x=solns(2);
plot(b0x,b0y,'g*','MarkerSize',10);
hold on
disp([b0x b0y]);

% Lets connect B0 with B1, B2, B3, and B4
B0B1x=[b0x,b1x];
B0B1y=[b0y,b1y];
plot(B0B1x,B0B1y,'r','LineWIdth',2);
hold on

B0B2x=[b0x,b2x];
B0B2y=[b0y,b2y];
plot(B0B2x,B0B2y,'y','LineWIdth',2);
hold on

B0B3x=[b0x,b3x];
B0B3y=[b0y,b3y];
plot(B0B3x,B0B3y,'m','LineWIdth',2);
hold on

B0B4x=[b0x,b4x];
B0B4y=[b0y,b4y];
plot(B0B4x,B0B4y,'g','LineWIdth',2);
hold on

A0B0x=[a0x,b0x];
A0B0y=[a0y,b0y];
plot(A0B0x,A0B0y,'k','LineWidth',2);
hold on

fprintf('***************************************GRAPHICAL RESULTS*****************************************\n');
A0A1=((a0x-a1x)^2+(a0y-a1y)^2)^0.5;
A1P1=((a1x-p1x)^2+(a1y-p1y)^2)^0.5;
A1B1=((a1x-b1x)^2+(a1y-b1y)^2)^0.5;
P1B1=((p1x-b1x)^2+(p1y-b1y)^2)^0.5;
B0B1=((b0x-b1x)^2+(b0y-b1y)^2)^0.5;
A0B0=((a0x-b0x)^2+(a0y-b0y)^2)^0.5;
fprintf('A0A1=%f\tA1P1=%f\tA1B1=%f\tP1B1=%f\tB0B1=%f\tA0B0=%f\n',A0A1,A1P1,A1B1,P1B1,B0B1,A0B0);

A0A2=((a0x-a2x)^2+(a0y-a2y)^2)^0.5;
A2P2=((a2x-p2x)^2+(a2y-p2y)^2)^0.5;
A2B2=((a2x-b2x)^2+(a2y-b2y)^2)^0.5;
P2B2=((p2x-b2x)^2+(p2y-b2y)^2)^0.5;
B0B2=((b0x-b2x)^2+(b0y-b2y)^2)^0.5;
A0B0=((a0x-b0x)^2+(a0y-b0y)^2)^0.5;
fprintf('A0A2=%f\tA2P2=%f\tA2B2=%f\tP2B2=%f\tB0B2=%f\tA0B0=%f\n',A0A2,A2P2,A2B2,P2B2,B0B2,A0B0);

A0A3=((a0x-a3x)^2+(a0y-a3y)^2)^0.5;
A3P3=((a3x-p3x)^2+(a3y-p3y)^2)^0.5;
A3B3=((a3x-b3x)^2+(a3y-b3y)^2)^0.5;
P3B3=((p3x-b3x)^2+(p3y-b3y)^2)^0.5;
B0B3=((b0x-b3x)^2+(b0y-b3y)^2)^0.5;
A0B0=((a0x-b0x)^2+(a0y-b0y)^2)^0.5;
fprintf('A0A3=%f\tA3P3=%f\tA3B3=%f\tP3B3=%f\tB0B3=%f\tA0B0=%f\n',A0A3,A3P3,A3B3,P3B3,B0B3,A0B0);

A0A4=((a0x-a4x)^2+(a0y-a4y)^2)^0.5;
A4P4=((a4x-p4x)^2+(a4y-p4y)^2)^0.5;
A4B4=((a4x-b4x)^2+(a4y-b4y)^2)^0.5;
P4B4=((p4x-b4x)^2+(p4y-b4y)^2)^0.5;
B0B4=((b0x-b4x)^2+(b0y-b4y)^2)^0.5;
A0B0=((a0x-b0x)^2+(a0y-b0y)^2)^0.5;
fprintf('A0A4=%f\tA4P4=%f\tA4B4=%f\tP4B4=%f\tB0B4=%f\tA0B0=%f\n',A0A4,A4P4,A4B4,P4B4,B0B4,A0B0);

% fill([a1x,b1x,p1x],[a1y,b1y,p1y],'r');
% hold on
% fill([a2x,b2x,p2x],[a2y,b2y,p2y],'y');
% hold on
% fill([a3x,b3x,p3x],[a3y,b3y,p3y],'m');
% hold on
% fill([a4x,b4x,p4x],[a4y,b4y,p4y],'g');
% hold on
axis('equal')