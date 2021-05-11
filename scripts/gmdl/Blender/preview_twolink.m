clr;
%% load data
load('deformData.mat');
Length0 = 65e-3;
FilterCoeff = 50;
Elong = smooth(1+epsilon,FilterCoeff);
Theta = smooth((angle-angle(1)),FilterCoeff);

%% preview
obj = Gmodel('SoftActuatorPlanarRedux.stl','TextureStretch',0.95);
obj = Blender(obj,'Scale',Length0);
obj.fix();
obj.bake.render();

f = gcf;
for ii = 1:fps(t,13):500
    
    obj.reset();
    obj = Blender(obj,'Curve',{'PCC+',Theta(ii),0,Elong(ii)});
    obj = Blender(obj,'Rotate',{'x',90});
    
    obj.update();
    
    axis tight;
    axis equal;
    axis([-0.05 0.05 -0.1 0]);
    %view(0,10);
    set(f,'Name',['TIME = ',num2str(t(ii),3),'(s)']);
end

