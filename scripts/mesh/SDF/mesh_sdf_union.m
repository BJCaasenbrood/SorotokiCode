clr;
%% generate SDF class
sdf1 = sRectangle(-0.5,0.5,-0.5,0.5);
sdf2 = sCircle(0, 0, 0.5);

sdf1 = sdf1.rotate(pi/4);
%% math operations
figure(101);
subplot(1,3,1);
S = sdf1;
sdf2.show();

subplot(1,3,2);
S = sdf2;
sdf1.show();

subplot(1,3,3);
S = sdf1.smoothunion(sdf2,0.2);
S.show()