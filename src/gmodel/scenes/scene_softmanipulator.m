function [base, link, grip] = scene_softmanipulator(Y1,Y2,varargin)
% RIG_SOFTMANIPULATOR set a predefined rig for a soft manipulator + gripper
L = 60;
Modes = [0,1,1,1,0,0];

if ~isempty(varargin)
   for ii = 1:2:numel(varargin) 
       if strcmp(varargin{2*ii-1},'L')
           L = varargin{2*ii};
       elseif strcmp(varargin{2*ii-1},'Modes')
           Modes = varargin{2*ii};
       end
   end
end

H1 = Gmodel('FrameSRMHolder.stl','Shading','Face'); 
H1 = Blender(H1,'Scale',(L/64));

H1 = Blender(H1,'Scale',{'z',-1});
L1 = Gmodel('Pneulink.stl');    
C0 = Gmodel('ConnectorRedux.stl','Shading','Face');
C1 = Gmodel('ConnectorRedux.stl','Shading','Face');
C2 = Gmodel('ConnectorRedux.stl','Shading','Face');
C2 = Blender(C2,'Scale',{'z',0.2});
C2 = Blender(C2,'Translate',{'z',0.125});
G1 = Gmodel('GripperFingerRedux.stl');

%SLAPRINT = (elasticresin);
SLAPRINT = (bluebase);

H1.Texture = egg;
L1.Texture = SLAPRINT;
G1.Texture = L1.Texture;
C1.Texture = H1.Texture;
C0.Texture = H1.Texture;

% shapes
shp1 = Shapes(Y1,Modes,'L0',L);
shp2 = Shapes(Y2,[0,0,size(Y2,2),0,0,0],'L0',L);

% complete rig
link = Rig(@(x) shp1.string(x),'Domain',L,'XTangent',true);
grip = Rig(@(x) shp2.string(x),'Domain',L);

link = link.add(L1,C1,C2);
link = link.parent(1,0,0);
link = link.parent(1,1,1);
link = link.parent(2,1,1);
link = link.parent(3,1,2);

link = link.texture(1,SLAPRINT);
link = link.texture(2,egg);
link = link.texture(3,SLAPRINT);

link.g0 = SE3(roty(pi),zeros(3,1));

base  = H1.render();
base0 = C0.render();

grip = grip.add(G1);
grip = grip.parent(1,0,0);
grip = grip.parent(1,1,1);
grip = grip.cpattern(1,{'z',3});
grip = grip.texture(1,SLAPRINT);

grip.g00 = SE3(rotx(deg2rad(-30)),[0,14,-7.5]);
link = link.computeFK(zeros(sum(Modes),1)); 
Ge = link.g(:,:,end);

grip.g0 = Ge(:,:,end);
grip    = grip.computeFK(zeros(size(Y2,2),1));

grip = grip.render();
link = link.render();

view(30,10);
axis tight;

end

