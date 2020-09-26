p1=[-2,8];
a1=[-4.724045618509668,7.739714239541714];
b1=[-4.837515403760070,10.862990065686912];
a0=[0.072253434046299,6.355493076991842];
b0=[-22.746103404754226,61.591109497526310];
plot(p1(1),p1(2),'r*', 'MarkerSize', 10);
hold on
plot(a1(1),a1(2),'r*', 'MarkerSize', 10);
hold on
plot(b1(1),b1(2),'r*', 'MarkerSize', 10);
hold on
plot(a0(1),a0(2),'r*', 'MarkerSize', 10);
hold on
plot(b0(1),b0(2),'r*', 'MarkerSize', 10);
hold on
A0A1x=[a0(1),a1(1)];
A0A1y=[a0(2),a1(2)];
plot(A0A1x,A0A1y,'m','Linewidth',2);
hold on
P1A1x=[p1(1),a1(1)];
P1A1y=[p1(2),a1(2)];
plot(P1A1x,P1A1y,'m','Linewidth',2);
hold on
B1A1x=[b1(1),a1(1)];
B1A1y=[b1(2),a1(2)];
plot(B1A1x,B1A1y,'m','Linewidth',2);
hold on
B1P1x=[b1(1),p1(1)];
B1P1y=[b1(2),p1(2)];
plot(B1P1x,B1P1y,'m','Linewidth',2);
hold on
B1B0x=[b1(1),b0(1)];
B1B0y=[b1(2),b0(2)];
plot(B1B0x,B1B0y,'m','Linewidth',2);
hold on
A0B0x=[a0(1),b0(1)];
A0B0y=[a0(2),b0(2)];
plot(A0B0x,A0B0y,'m','Linewidth',2);
hold on