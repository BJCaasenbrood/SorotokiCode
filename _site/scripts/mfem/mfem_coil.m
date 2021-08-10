clr;
%% set signed distance function
D = 5.0;
R = 4.0;
H = 5.0;
W = 10;

sdf =  cUnion(cHelix(D,0,0,R,H,W),...
       cUnion(cHelix(D*cos(1.3333*pi),D*sin(1.3333*pi),0,R,H,W),...
              cHelix(D*cos(0.6666*pi),D*sin(0.6666*pi),0,R,H,W)));

%% generate mesh
mag = Mmesh(sdf,'NElem',1000);
mag = mag.generate();

%% compute inductance
mag = Blender(mag,'Curve',{'PCC',20,0.0});
mag.inductance();

%% render
mag.render();

