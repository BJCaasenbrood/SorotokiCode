clc; clear; close all;
%% parameters
H = 25; 
W = 3; 

Dist = @(X) dCube(X,-W/2,W/2,-W/2,W/2,0,H);
obj = Gmodel(Dist,[-W/2,W/2,-W/2,W/2,.1,H-.1]);

%% set texture
obj.Texture = base;
obj.render(); 

%% blend shapes
Frame = 75;

v0 = .5;
v1 = .5;
a0 = 0;
a1 = 0;

for ii = 0:1
obj = obj.resetNode();
deg = polyq5((ii/Frame),0,1,v0,v1,a0,a1);
obj = Blender(obj,'Rotate',{'z',360*deg});
obj = Blender(obj,'Twist',{'z',360*1.25});
obj = Blender(obj,'Curve',{'PCC+',361.5,0,1.15});
axis tight;
view(10,10);
axis([-1.5 10.5 -1.5 1.5 -6 6]);
obj.update();
if ii == 1, gif('bla.gif','frame',gcf,'DelayTime',1/24); 
elseif ii > 1, gif; end
fprintf('frame number: %1.0f \n',ii);
end


