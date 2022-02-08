clr;
%% loading duck
obj = Gmodel('rubberducky.stl');
obj = Blender(obj,'Scale',1/160);
obj = Blender(obj,'Translate',{'z',-0.25});
obj = Blender(obj,'Rotate',{'y',15});
obj.Texture = diffuse(0.15); 
obj.bake.render();

%% setting up soft kinematics
N = 150;
M = 5;
x = linspace(0,1,N).';
Y = zeros(N,M);

for ii = 1:M
    Y(:,ii) = chebyshev(x,ii-1); % legendre 
end

%Y = gsogpoly(Y);
shp = Shapes(Y,[0,M,0,0,0,0]);

%% setting up soft arm
gmdl = Gmodel('Arm.stl');
gmdl = gmdl.set('Emission', [0.9 0.8 0.8],...
    'SSSPower',0.1,'SSSRadius',0.15,'SSS',true);

rig = Rig(@(x) shp.string(x),...
    'Domain',1);

rig = rig.add(gmdl);
rig = rig.parent(1,0,0);
rig = rig.parent(1,1,.999);
rig = rig.render();

rig.g0 = [rot2quat(rotx(pi/2)),-0.2,0.3,0];

rig = rig.computeFK(zeros(shp.NDim,1));

axis tight;
view(30,30);
obj.update();
rig.update();

%% set grab point
Sd = sCircle(0.3,0.3,0.05);

%% optimization
q = 1e-3*rand(shp.NDim,1);
% qopt = fminunc(@(x) Objective(x,shp,Sd),q);
% 
% rig = rig.computeFK(qopt(:));
% rig = rig.update();
e = 0.1;

%% functions
while norm(e) > 1e-4
     
     [g,J] = shp.string(q);
    
     p = reshape(g(1:3,4,:),3,[]).';
     
     [XY,~]  = ClosestPointOnSDF(Sd,p(:,[1 3]));
     [T,N,B] = Sd.normal(XY);
     
     F = 0;  
     
     dq = 0;

     for ii = 1:shp.NNode
         
         pd = [XY(ii,1);0;XY(ii,2)];
         Rd = [T(ii,1),B(ii,1),N(ii,1);
               T(ii,3),B(ii,3),N(ii,3);
               T(ii,2),B(ii,2),N(ii,2)];
         
         gd = SE3(flipud(Rd),pd);
         
         dF = EnergyController(g(:,:,ii),gd);
         
         dq = dq + 1e-2*J(:,:,ii).'*dF;
         %F = F + shp.Sigma(ii)*(dF).'*dF;
     end
     
     q = q + dq;
     e = dF.'*dF;
    
     rig = rig.computeFK(q(:));
     rig = rig.update();
     
     %Cost = F;
end

function Cost = Objective(q,shp,Sd)
    
     [g,~] = shp.string(q);
    
     p = reshape(g(1:3,4,:),3,[]).';
     
     [XY,~]  = ClosestPointOnSDF(Sd,p(:,[1 3]));
     [T,N,B] = Sd.normal(XY);
     
     F = 0;     

     for ii = 1:shp.NNode
         
         pd = [XY(ii,1);0;XY(ii,2)];
         Rd = [T(ii,1),B(ii,1),N(ii,1);
               T(ii,3),B(ii,3),N(ii,3);
               T(ii,2),B(ii,2),N(ii,2)];
         
         gd = SE3(flipud(Rd),pd);
         
         dF = EnergyController(g(:,:,ii),gd);
         F = F + shp.Sigma(ii)*(dF).'*dF;
     end
    
   Cost = F;
     
end
function Fu = EnergyController(g,gd)
    
    k1 = 0.05;
    k2 = 1;
    
    Kp = diag([k1,k1,k1,k2,k2,k2]);

    Xi = logmapSE3(g\gd);
    Fu = Kp*tmapSE3(Xi)*isomse3(Xi);
end
function [XY, D] = ClosestPointOnSDF(sdf,P)
    [XY, D] = distance2curve(sdf.Node,P);
end
