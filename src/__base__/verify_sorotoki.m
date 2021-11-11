function out = verify_sorotoki
%VERIFY_SOROTOKI Summary of this function goes here
%   Detailed explanation goes here
clr;

global dt;
dt = 0.01;
fprintf('* Starting verification check - SOROTOKI \n'); pause(0.2)
fprintf('* The verification check might take a couple of minutes. \n');
fprintf('  The analysis will perform a full check shortly... \n');
pause(1);

List = cell(0);

% check plotting within SOROTOKI
List = verifySorotoki(List,@(x) checkPlot,'Plotting tools'); 
close(101);

List = verifySorotoki(List,@(x) checkSDF,'Class Sdf.m');
close(101);

List = verifySorotoki(List,@(x) checkMesh,'Class Mesh.m');
close(101);

List = verifySorotoki(List,@(x) checkFem,'Class Fem.m');
close(101);


out = List;
end

function checkPlot
global dt
figure(101);
set(gcf,'Name','Checking Plotting tools');
subplot(2,2,1);
cout('\t col()\t');
bic;
colorshow();
boc;

subplot(2,2,2);
[X,Y,Z] = peaks(40);
surf(X,Y,Z,'linestyle','none')
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
imshow('FEMbot.png');
boc;
end

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

function checkMesh
global dt

F = [1,2,3,4];
V = [0,0;5,0;5,5;0,5];

sdf = @(x) dCircle(x,0,0,1);

cout('\t Mesh(V,F)'); bic; 
%msh1 = Mesh(V,F);
%msh1 = msh1.generate();
boc;

cout('\t Mesh(Sdf)'); bic;
msh2 = Mesh(sdf,'BdBox',[-1,1,-1,1],'NElem',250);
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

function checkFem
global dt

F = [1,2,3,4];
V = [0,0;5,0;5,5;0,5];

cout('\t Fem(Mesh)'); bic; 
msh = Mesh(V,F);
msh = msh.generate();
fem = Fem(msh,'TimeStep',1/5,'PrescribedDisplacement',true,'ShowProcess',false);
boc;

cout('\t Fem.AddConstraint'); bic(1); 
fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,0]);
fem = fem.AddConstraint('Support',fem.FindNodes('SW'),[1,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('Right'),[5,0]);
boc;

cout('\t Fem.Material'); bic; 
fem.Material = Ecoflex0050(0.1);
boc;

cout('\t Fem.solve()'); bic; 
[~,U] = fem.solve();
boc;


end

function list = verifySorotoki(list,exec,label)
list = verifyFunction(list,label,0);
try
    exec(0);
    list = verifyFunction(list,label,1);
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
    cout('red', input.message);
    cout('\n');
end

end

function bic(n)
if nargin < 1, cout('\t\t ');
else,  cout('\t ');
end
cout('key','Checking...'); 
end
function boc, cout('\b\b\b\b\b\b\b\b\b\b\b\b '); cout('green','Passed \n'); end