
% VERIFY_SOROTOKI is the installation checker, and ensures all the functions
% and classes of SOROTOKI are working accordingly.
%
% Usage:
%   x = verifySorotoki();   % calls the installation checker
%
%   x is a cell with Exception Catches -- showing the errrors that might 
%       occuring during the verification.
function vargout = verifySorotoki(input)

if nargin > 0
    clc;
else
    input = [];
end

global dt;
dt = 0.01;
fprintf('* Starting verification check - SOROTOKI \n'); pause(0.2)
fprintf('* The verification check might take a couple of minutes. \n');
fprintf('  The analysis will perform a full check shortly... \n');
pause(1);

List = cell(0);
    

if isa(input,'char')
    switch lower(input)
        case {'base','plot','plotting'}
             List = verifySorotokiSoftware(List,@(x) checkPlot,'Plotting tools');        
        case {'model','mdl'}
            List = verifySorotokiSoftware(List,@(x) checkModel,'Class Model.m');
        case {'sdf', 'signed distance'}
            List = verifySorotokiSoftware(List,@(x) checkSDF,'Class Sdf.m');
        case {'fem'}
            List = verifySorotokiSoftware(List,@(x) checkFem,'Class Fem.m');
        case {'mesh','msh'}
            List = verifySorotokiSoftware(List,@(x) checkMesh,'Class Mesh.m');
        case {'graphics','gmodel','graphic'}
            List = verifySorotokiSoftware(List,@(x) checkGmodel,'Class Gmodel.m');
        otherwise
            cout('Unknown src files\n');
    end
else  

    %check plotting within SOROTOKI
    List = verifySorotokiSoftware(List,@(x) checkPlot,'Plotting tools');
    
    List = verifySorotokiSoftware(List,@(x) checkSDF,'Class Sdf.m');
    
    List = verifySorotokiSoftware(List,@(x) checkMesh,'Class Mesh.m');
    
    List = verifySorotokiSoftware(List,@(x) checkFem,'Class Fem.m');
    
    List = verifySorotokiSoftware(List,@(x) checkGmodel,'Class Gmodel.m');
    
    List = verifySorotokiSoftware(List,@(x) checkModel,'Class Model.m');
    
end

vargout = List;
end
%---------------------------------------------------------- check Plotting
function checkPlot
global dt
figure(101);
set(gcf,'Name','Checking Plotting tools');
subplot(2,2,1);
cout('\t color\t');
bic;
colorshow();
boc;

subplot(2,2,2);
colorwheel();
% [X,Y,Z] = peaks(40);
% surf(X,Y,Z,'linestyle','none')
cout('\t colormap')
bic;
colormap(turbo); drawnow; pause(dt);
colormap(barney); drawnow; pause(dt);
colormap(blackwhite); drawnow; pause(dt);
colormap(bluesea); drawnow; pause(dt);
colormap(bounce); drawnow; pause(dt);
colormap(inferno); drawnow; pause(dt);
colormap(metro); drawnow; pause(dt);
colormap(noir); drawnow; pause(dt);
colormap(turbo); drawnow; pause(dt);
colormap(viridis); drawnow; pause(dt); 
boc;

subplot(2,2,3);
cout('\t cmapping()')
bic;
showColormap(turbo); drawnow; pause(dt);
showColormap(turbo(-1)); drawnow; pause(dt);
showColormap(turbo(0)); drawnow; pause(dt);
showColormap(turbo(-100)); drawnow; pause(dt);
showColormap(turbo(100)); drawnow; pause(dt);
boc;

subplot(2,2,4);
cout('\t verify images')
bic;
imshow('Pneunet.png');
boc;
end
%---------------------------------------------------------- check Sdf class
function checkSDF
global dt
f1 = @(x) sqrt((x(:,1)).^2+(x(:,2)).^2)-0.75;
f2 = @(x) sqrt((x(:,1)).^2+(x(:,2)).^2 + +(x(:,3)).^2)-0.75;
f3 = @(x) sqrt((x(:,1)+0.5).^2+(x(:,2)).^2 + +(x(:,3)).^2)-0.75;

cout('\t Sdf(...) 2D')
bic; sdf1 = Sdf(f1,'BdBox',[-1,1,-1,1]); boc;

cout('\t Sdf(...) 3D')
bic; 
sdf2 = Sdf(f2,'BdBox',[-1,1,-1,1,-1,1]); 
sdf3 = Sdf(f3,'BdBox',[-1.5,1.5,-1,1,-1,1]);
boc;

figure(101);
set(gcf,'Name','Checking Signed Distance Functions Sdf.m');
subplot(2,2,1);
cout('\t Sdf.show() 2D')
bic; sdf1.show(); boc;  drawnow; pause(2*dt);
subplot(2,2,2);
cout('\t Sdf.show() 3D')
bic; sdf2.show(); boc; drawnow; pause(2*dt);

subplot(2,2,3);
cout('\t sdf1 + sdf2 ');
bic; 
s = sdf2 + sdf3; 
s.show(); boc; 
drawnow; pause(5*dt);

cout('\t sdf1 - sdf2 ');
bic; 
s = sdf2 - sdf3; cla;
s.show(); boc; 
drawnow; pause(5*dt);

subplot(2,2,4);
cout('\t sdf1/sdf2 ');
bic; 
s = sdf2/sdf3; cla;
s.show(); boc; 
drawnow; pause(5*dt);
end
%--------------------------------------------------------- check Mesh class
function checkMesh
global dt

F = [1,2,3,4];
V = [0,0;5,0;5,5;0,5];

sdf = sCircle(1);

cout('\t Mesh(V,F)'); bic; 
msh1 = Mesh(V,F);
msh1 = msh1.generate();
boc;

cout('\t Mesh(Sdf)'); bic;
msh2 = Mesh(sdf,'BdBox',[-1,1,-1,1],'NElem',50);
msh2 = msh2.generate(); 
boc;

cout('\t Mesh.show()'); bic;
figure(101);
msh1.show();
drawnow; pause(5*dt);
msh2.show();
drawnow; pause(5*dt);
boc;
end
%---------------------------------------------------------- check FEM class
function checkFem
sdf = sRectangle(0,2,0,1);

cout('\t Fem(Mesh)'); bic; 
msh = Mesh(sdf,'Quads',[8,3]);
msh = msh.generate();
fem = Fem(msh,'ShowProcess',false);
boc;

cout('\t Fem.show()'); bic; 
boc;

cout('\t Fem.AddConstraint()'); bic(1); 
fem = fem.addSupport(fem.FindNodes('Left'),[1,1]);
fem = fem.addSupport(fem.FindNodes('Right'),[0,1]);
fem = fem.addDisplace(fem.FindNodes('Right'),[1,0]);
boc;

fem.Material = Ecoflex0050(10);

cout('\t Fem.solve()'); bic; 
fem.solve();
boc;

cout('\t NeoHookean()'); bic; 
fem = fem.reset();
fem.Material = NeoHookeanMaterial(); fem.solve();
boc();

cout('\t MooneyMaterial()'); bic(1); 
fem = fem.reset();
fem.Material = MooneyMaterial(); fem.solve();
boc();

cout('\t Yeoh() '); bic; 
fem = fem.reset();
fem.Material = YeohMaterial(); fem.solve();
boc();

cout('\t Fem.optimize()'); bic; 
sdf = sRectangle(0,5,0,2);
msh = Mesh(sdf,'NElem',600);
msh = msh.generate();
fem = Fem(msh,'ShowProcess',0,'Nonlinear',0,...
    'ChangeMax',0.2,'MaxIterationMMA',30);

fem = fem.set('FilterRadius',0.1,'VolumeInfill',0.2,...
              'Penal',3,'ReflectionPlane',[-1,0]);
        
fem = fem.addSupport(fem.FindNodes('SE'),[0,1]);
fem = fem.addSupport(fem.FindNodes('Left'),[1,0]);
fem = fem.addLoad(fem.FindNodes('Location',[0,2],1),[0,-1e-3]);

fem = fem.initialTopology('Hole',[1,1;3,1],0.5);
fem.Material = NeoHookeanMaterial(1,0.33);
fem.optimize();
fem.show('ISO',0.25);
boc;

end
%---------------------------------------------------------- check FEM class
function checkGmodel
global dt
cout('\t Gmodel(.stl)'); bic; 
obj = Gmodel('Bunny.stl','ShowProcess',false); 
boc;

cout('\t Gmodel.bake()'); bic; 
obj = obj.bake(); boc;

cout('\t Gmodel.render()'); bic(1); 
obj = obj.render(); view(20,20); drawnow; pause(2*dt);
boc

cout('\t Gmodel.update()'); bic(1); 
shading = linspace(0,1,10);
for ii = shading
    obj.Texture = diffuse(ii);
    obj.update(); 
end
boc;

cout('\t Scale \t'); bic(); 
    obj1 = Blender(obj,'Scale',{'z',1.5});
    obj1.update; drawnow;
boc;     
    
cout('\t Rotate '); bic(); 
    obj1 = Blender(obj,'Rotate',{'y',50});
    obj1.update; drawnow;
boc;    

cout('\t Twist \t'); bic(); 
    obj = obj.reset;
    obj1 = Blender(obj,'Twist',{'z',220});
    obj1.update; drawnow;
boc; 
    

end
%---------------------------------------------------------- check FEM class
function checkModel
global dt
M = 5;
N = 200;

X = linspace(0,1,N).';
Y = [];

cout('\t pcc()  '); bic; 
for ii = 1:M
    Y(:,ii) = pcc(X,ii,M);
end
boc;
cout('\t pwl()  '); bic; 
for ii = 1:M
    Y(:,ii) = pwl(X,ii,M);
end
boc;

cout('\t legendre()'); bic; 
for ii = 1:M
    Y(:,ii) = legendre(X,ii-1); 
end
boc;

cout('\t chebyshev()'); bic; 
for ii = 1:M
    Y(:,ii) = chebyshev(X,ii-1); 
end
boc;

cout('\t gsogpoly()'); bic; 
Y = gsogpoly(Y,X);
boc;

cout('\t Shapes(Y,Modes)'); bic(1); 
shp = Shapes(Y,[0,M,0,0,0,0],'L0',100);   
boc;

cout('\t Shapes.string()'); bic(1); 
figure(101);
Q = 2e-2*diag(1:M);

for ii = 1:M
   cla;
   g = shp.string(Q(:,ii));
   plotSE2(g,'xz',0.1); axis equal; 
   axis(100*[0,1,-1,1]); axis off;
   axis tight;
   drawnow;
   pause(5*dt);
end
boc();

cout('\t Shapes.string()'); bic(1); 
p = shp.FK(Q(:,ii));
plot(p(:,1),p(:,3),'LineW',3);
boc();

E    = 3.00;     % Young's modulus in Mpa
Nu   = 0.42;     % Poisson ratio
shp.Material = NeoHookeanMaterial(E,Nu);
shp.Material.Zeta = 0.03;
shp.Gravity = [9.81e3; 0; 0];

cout('\t Shapes.rebuild()'); bic(1); 
shp = shp.rebuild(); boc();

cout('\t Model(shp) \t'); bic(1); 
mdl = Model(shp,'TimeStep',1/120,'TimeEnd',1,'ShowProcess',false);
mdl.X0(1) = 0.12;
boc();

cout('\t Model.simulate()'); bic(1); 
mdl = mdl.simulate(); 
boc;

cout('\t Rig(Shape) \t'); bic(1); 
gmdl = Gmodel('Arm.stl','ShowProcess',false);
rig    = Rig(@(x) shp.string(x),'Domain',100);
rig    = rig.add(gmdl);
rig    = rig.parent(1,0,0);
rig    = rig.parent(1,1,1);
rig    = rig.texture(1,softmath);
rig.g0 = SE3(roty(-pi),zeros(3,1));

figure(101); clf;
rig = rig.render();
axis(100*[-.5 .5 -.5 .5 -1 0.1]);
view(30,30);
rig = rig.computeFK(mdl.Log.x(ii,1:M));
rig = rig.update();

for ii = 1:fps(mdl.Log.t,70):length(mdl.Log.x)

    rig = rig.computeFK(mdl.Log.x(ii,1:M));
    rig = rig.update();
    
    axis(100*[-.5 .5 -.5 .5 -1 0.1]);
    view(30,30);
    drawnow();
end

boc();

end
%--------------------------------------------------------------------------
function list = verifySorotokiSoftware(list,exec,label)
list = verifyFunction(list,label,0);
try
    exec(0);
    list = verifyFunction(list,label,1);
    close(101);
catch e
    list = verifyFunction(list,e,-1);
end

end
function List = verifyFunction(List,input,index)
if index == 0
    cout('text','* ');
    cout('text','Assesing: ');
    str = char(input);
    cout('key',str);
    cout('key','... \n');
    List{end+1,1} = str;
elseif index == 1
    List{end,2} = [];
    cout(' ');
    cout('green','\t Complete!\n');
else
    List{end,2} = input;
    cout(' ');
    cout('red','\t Unsuccesfull!\n');
    %cout('red', input.message);
    cout('\n');
end

end
function bic(n)
if nargin < 1
       cout('\t\t ');
else,  cout('\t ');
end
cout('key','Checking...'); 
end
function boc
cout('\b\b\b\b\b\b\b\b\b\b\b\b '); 
cout('green','Passed \n'); 
end