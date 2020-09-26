% Analytical technique for path generation problem for three position with three prescribed timing

% Points through which a point on coupler will pass
p1=1+1i;%p1=input('Eneter first point through which a point on coupler point will pass=');
p2=-0.4+0.24i;%p2=input('Eneter second point through which a point on coupler point will pass=');
p3=0-1.3i;%p3=input('Eneter third point through which a point on coupler point will pass=');

% Lets find DYAD for second and third position with respect to first position
dyad2=p2-p1;
dyad3=p3-p1;
disp(dyad2);
disp(dyad3);

% Rotational angle of crank with respect to different position
thetal12=126;%theta12=input('Enter crank rotation angle in degrees when coupler moves from point P1 to point P2');
theta12=thetal12*pi/180; %%degree to radian angle
A=exp(theta12*1i)-1;
disp(A);

thetal13=252;%theta13=input('Enter crank rotation angle in degrees when coupler moves from point P1 to point P3');
theta13=thetal13*pi/180; %% degree to radian angle
B=exp(theta13*1i)-1;
disp(B);

% Assumed angle between the input and coupler link for different postion
alphal12=-6;%alphal12=input('Enter angle between input link and coupler link in degrees when coupler moves from point P1 to point P2');
alpha12=alphal12*pi/180; %%degree to radian angle 
C=exp(alpha12*1i)-1;
disp(C);

alphal13=37;%alphal13=input('Enter angle between input link and coupler link in degrees when coupler moves from point P1 to point P3');
alpha13=alphal13*pi/180; %% degree to radian angle
D=exp(alpha13*1i)-1;
disp(D);

% Assumed angle between the output and coupler link for different postion
fil12=33;%fil12=input('Enter angle between output link and coupler link in degrees when coupler moves from point P1 to point P2');
fi12=fil12*pi/180; %%degree to radian angle
E=exp(fi12*1i)-1;
disp(E);

fil13=37;%fil13=input('Enter angle between output link and coupler link in degrees when coupler moves from point P1 to point P3');
fi13=fil13*pi/180; %% degree to radian angle
F=exp(fi13*1i)-1;
disp(F);


% Lets find link length and its orientation with horizontal
z1=det([dyad2 C;dyad3 D])/det([A C;B D]);
disp(z1);
z2=det([A dyad2;B dyad3])/det([A C;B D]);
disp(z2);
z3=det([dyad2 C;dyad3 D])/det([E C;F D]);
disp(z3);
z4=det([E dyad2;F dyad3])/det([E C;F D]);
disp(z4);
z5=z2-z4;
disp(z5);
z6=z1+z5-z3;
disp(z6);

% Lets convert the link vectors into polar form
[theta1,r1]=cart2pol(real(z1),imag(z1));
disp([theta1,r1]);
[theta2,r2]=cart2pol(real(z2),imag(z2));
disp([theta2,r2]);
[theta3,r3]=cart2pol(real(z3),imag(z3));
disp([theta3,r3]);
[theta4,r4]=cart2pol(real(z4),imag(z4));
disp([theta4,r4]);
[theta5,r5]=cart2pol(real(z5),imag(z5));
disp([theta5,r5]);
[theta6,r6]=cart2pol(real(z6),imag(z6));
disp([theta6,r6]);

% Lets plot the mechanism
a0x=0;
a0y=0;
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