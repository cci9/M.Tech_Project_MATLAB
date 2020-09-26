 %%%%%%%%% Path Generation Three Precision Point Problem with Prescribed Timing %%%%%%%%%%
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

% Rotational angle of crank with respect to different position
thetal12=60;%theta12=input('Enter crank rotation angle in degrees when coupler moves from point P1 to point P2');
thetal21=thetal12;
theta21=thetal21*pi/180; %%degree to radian angle 

thetal13=110;%theta13=input('Enter crank rotation angle in degrees when coupler moves from point P1 to point P3');
thetal31=thetal13;
theta31=thetal31*pi/180; %% degree to radian angle 

% Connect different links
A0B0x=[a0x,b0x];
A0B0y=[a0y,b0y];
A0P1x=[a0x,p1x];
A0P1y=[a0y,p1y];
A0P2x=[a0x,p2x];
A0P2y=[a0y,p2y];
A0P3x=[a0x,p3x];
A0P3y=[a0y,p3y];

% Mechanism Ploting
plot(A0B0x,A0B0y,'green','LineWidth',2);
hold on
plot(A0P1x,A0P1y,'--');
hold on
plot(A0P2x,A0P2y,'--');
hold on
plot(A0P3x,A0P3y,'--');
hold on

% Find the length and angle of each position point with fixed point
L1=((p1x-a0x)^2+(p1y-a0y)^2)^0.5;
alphaa0p1=atan((p1y-a0y)/(p1x-a0x));

L2=((p2x-a0x)^2+(p2y-a0y)^2)^0.5;
alphaa0p2=atan((p2y-a0y)/(p2x-a0x));

L3=((p3x-a0x)^2+(p3y-a0y)^2)^0.5;
alphaa0p3=atan((p3y-a0y)/(p3x-a0x));

% Find the different poition of P2, P3 after rotation by thetal21 and thetal31
p2nx=a0x+L2*cos(theta21+alphaa0p2); %p2xn represents new posiotion of P2 point
p2ny=a0y+L2*sin(theta21+alphaa0p2);
A0P2nx=[a0x,p2nx];
A0P2ny=[a0y,p2ny];
plot(A0P2nx,A0P2ny,'--');
hold on
L2new=((p2nx-a0x)^2+(p2ny-a0y)^2)^0.5;

p3nx=a0x+L3*cos(theta31+alphaa0p3); %p2xn represents new posiotion of P2 point
p3ny=a0y+L3*sin(theta31+alphaa0p3);
A0P3nx=[a0x,p3nx];
A0P3ny=[a0y,p3ny];
plot(A0P3nx,A0P3ny,'--');
hold on
L3new=((p3nx-a0x)^2+(p3ny-a0y)^2)^0.5;

% Connect new position of P2 with P1 and new poition of P3 with new position of P2
P1xP2nx=[p1x,p2nx];
P1yP2ny=[p1y,p2ny];
plot(P1xP2nx,P1yP2ny,'--');
hold on
P2nxP3nx=[p2nx,p3nx];
P2nyP3ny=[p2ny,p3ny];
plot(P2nxP3nx,P2nyP3ny,'--');
hold on

% Finding perpendicular bisector of P2nP1 and P2nP3n
p1p2nxCenter=(p1x+p2nx)/2;
p1p2nyCenter=(p1y+p2ny)/2;
plot(p1p2nxCenter, p1p2nyCenter, 'r*', 'MarkerSize', 3);
hold on

p2np3nxCenter=(p2nx+p3nx)/2;
p2np3nyCenter=(p2ny+p3ny)/2;
disp(p2np3nxCenter);
disp(p2np3nyCenter);
plot(p2np3nxCenter, p2np3nyCenter, 'r*', 'MarkerSize', 3);
hold on

mp1p2n=(p1y-p2ny)/(p1x-p2nx);
mp1p2n_ortho=-(1/mp1p2n);
y_ortho1=mp1p2n_ortho*(x-p1p2nxCenter)+p1p2nyCenter;

mp2np3n=(p2ny-p3ny)/(p2nx-p3nx);
mp2np3n_ortho=-(1/mp2np3n);
y_ortho2=mp2np3n_ortho*(x-p2np3nxCenter)+p2np3nyCenter;
dx=0.1;
x=-20:dx:20;
plot(x,y_ortho1);
hold on
plot(x,y_ortho2);
hold on

% Intersection point
y_ortho1=mp1p2n_ortho*(x-p1p2nxCenter)+p1p2nyCenter;
y_ortho2=mp2np3n_ortho*(x-p2np3nxCenter)+p2np3nyCenter;
eqns_function=@(var)[var(1)-(mp1p2n_ortho)*(var(2)-p1p2nxCenter)-p1p2nyCenter;
                     var(1)-(mp2np3n_ortho)*(var(2)-p2np3nxCenter)-p2np3nyCenter];
initials=[0;0];
solns=fsolve(eqns_function,initials,...
             optimoptions('fsolve','Display','off'))       
a1y=solns(1);
a1x=solns(2);
plot(a1x,a1y,'r*','MarkerSize',3);
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
fprintf('\n\t\t\t\t\t\t\t\t****GRAPHICAL RESULTS****\n');
A0A1Length=((a0x-a1x)^2+(a0y-a1y)^2)^0.5;
A1P1Length=((a1x-p1x)^2+(a1y-p1y)^2)^0.5;
P1B1Length=((p1x-b1x)^2+(p1y-b1y)^2)^0.5;
A1B1Length=((a1x-b1x)^2+(a1y-b1y)^2)^0.5;
B1B0Length=((b1x-b0x)^2+(b1y-b0y)^2)^0.5;
A0B0Length=((a0x-b0x)^2+(a0y-b0y)^2)^0.5;
fprintf('A0A1=%f\tA1P1=%f\tP1B1=%f\tA1B1=%f\tB1B0=%f\tA0B0=%f\n',A0A1Length,A1P1Length,P1B1Length,A1B1Length,B1B0Length,A0B0Length);

A0A2Length=((a0x-a2x)^2+(a0y-a2y)^2)^0.5;
A2P2Length=((a2x-p2x)^2+(a2y-p2y)^2)^0.5;
P2B2Length=((p2x-b2x)^2+(p2y-b2y)^2)^0.5;
A2B2Length=((a2x-b2x)^2+(a2y-b2y)^2)^0.5;
B2B0Length=((b2x-b0x)^2+(b2y-b0y)^2)^0.5;
fprintf('A0A2=%f\tA2P2=%f\tP2B2=%f\tA2B2=%f\tB2B0=%f\tA0B0=%f\n',A0A2Length,A2P2Length,P2B2Length,A2B2Length,B2B0Length,A0B0Length);

A0A3Length=((a0x-a3x)^2+(a0y-a3y)^2)^0.5;
A3P3Length=((a3x-p3x)^2+(a3y-p3y)^2)^0.5;
P3B3Length=((p3x-b3x)^2+(p3y-b3y)^2)^0.5;
A3B3Length=((a3x-b3x)^2+(a3y-b3y)^2)^0.5;
B3B0Length=((b3x-b0x)^2+(b3y-b0y)^2)^0.5;
fprintf('A0A3=%f\tA3P3=%f\tP3B3=%f\tA3B3=%f\tB3B0=%f\tA0B0=%f\n\n',A0A3Length,A3P3Length,P3B3Length,A3B3Length,B3B0Length,A0B0Length);

axis([-10 10 0 12])
axis('equal')

alpha1=atan((p1y-a1y)/(p1x-a1x));
if alpha1<0
    alpha1=pi+alpha1;
end

alpha2=atan((p2y-a2y)/(p2x-a2x));
if alpha2<0
    alpha2=pi+alpha2;
end

alpha3=atan((p3y-a3y)/(p3x-a3x));
if alpha3<0
    alpha3=alpha3+pi;
end

fi1=atan((b1y-b0y)/(b1x-b0x));
if fi1<0
    fi1=pi+fi1;
end

fi2=atan((b2y-b0y)/(b2x-b0x));
if fi2<0
    fi2=fi2+pi;
end

fi3=atan((b3y-b0y)/(b3x-b0x));
if fi3<0
    fi3=fi3+pi;
end

% Points through which a point on coupler will pass
p1=3+8i;%p1=input('Eneter first point through which a point on coupler point will pass=');
p2=4+8.5i;%p2=input('Eneter second point through which a point on coupler point will pass=');
p3=5+8.2i;%p3=input('Eneter third point through which a point on coupler point will pass=');

% Lets find DYAD for second and third position with respect to first position
dyad2=p2-p1;
dyad3=p3-p1;

% Rotational angle of crank with respect to different position
thetal12=-thetal12;%theta12=input('Enter crank rotation angle in degrees when coupler moves from point P1 to point P2');
theta12=thetal12*pi/180; %%degree to radian angle
A=exp(theta12*1i)-1;

thetal13=-thetal13;%theta13=input('Enter crank rotation angle in degrees when coupler moves from point P1 to point P3');
theta13=thetal13*pi/180; %% degree to radian angle
B=exp(theta13*1i)-1;

% Assumed angle between the input and coupler link for different postion
alpha12=alpha2-alpha1; %%degree to radian angle 
C=exp(alpha12*1i)-1;

alpha13=alpha3-alpha1; %% degree to radian angle
D=exp(alpha13*1i)-1;

% Assumed angle between the output and coupler link for different postion
fi12=fi2-fi1; %%degree to radian angle
E=exp(fi12*1i)-1;

fi13=fi3-fi1; %% degree to radian angle
F=exp(fi13*1i)-1;

% Lets find link length and its orientation with horizontal
z1=det([dyad2 C;dyad3 D])/det([A C;B D]);
% disp(z1);
z2=det([A dyad2;B dyad3])/det([A C;B D]);
% disp(z2);
z3=det([dyad2 C;dyad3 D])/det([E C;F D]);
% disp(z3);
z4=det([E dyad2;F dyad3])/det([E C;F D]);
% disp(z4);
z5=z2-z4;
% disp(z5);
z6=z1+z5-z3;
% disp(z6);

% Lets convert the link vectors into polar form
[theta1,r1]=cart2pol(real(z1),imag(z1));
% disp([theta1,r1]);
[theta2,r2]=cart2pol(real(z2),imag(z2));
% disp([theta2,r2]);
[theta3,r3]=cart2pol(real(z3),imag(z3));
% disp([theta3,r3]);
[theta4,r4]=cart2pol(real(z4),imag(z4));
% disp([theta4,r4]);
[theta5,r5]=cart2pol(real(z5),imag(z5));
% disp([theta5,r5]);
[theta6,r6]=cart2pol(real(z6),imag(z6));
% disp([theta6,r6]);

% Lets dispay analytical result on the commant window
fprintf('\t\t\t\t\t\t\t\t****ANALYTICAL RESULTS****\n');
fprintf('A0A1=%f\tA1P1=%f\tP1B1=%f\tA1B1=%f\tB1B0=%f\tA0B0=%f\n',r1,r2,r4,r5,r3,r6);

% Lets plot the mechanism
a0x=-2;
a0y=2;
a1x=a0x+(r1*cos(theta1));
a1y=a0y+(r1*sin(theta1));
plot([a0x a1x],[a0y a1y])
hold on
p1x=a1x+(r2*cos(theta2));
p1y=a1y+(r2*sin(theta2));
plot([a1x p1x],[a1y p1y])
hold on
b0x=a0x+(r6*cos(theta6));
b0y=a0y+(r6*sin(theta6));
plot([a0x b0x],[a0y b0y])
hold on
b1x=b0x+(r3*cos(theta3));
b1y=b0y+(r3*sin(theta3));
plot([b0x b1x],[b0y b1y])
hold on
p1nx=b1x+(r4*cos(theta4));
p1ny=b1y+(r4*sin(theta4));
plot([b1x p1nx],[b1y p1ny])
hold on
b1nx=a1x+(r5*cos(theta5));
b1ny=a1y+(r5*sin(theta5));
plot([a1x b1nx],[a1y b1ny])
hold on