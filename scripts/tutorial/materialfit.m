% TUTORIAL 1: Material fit demo
clr;
Data = load('unknownMaterialSamples.mat');
%% assigning dimenions and elongation
N  = 30;
H  = 20;  % height of specimen
W  = 20;  % width of specimen
dL = H*5; % elongation of specimen
 
%% generate data set
e   = linspace(0,dL/H,N);
svm = exp(0.5*e)-1;

%% signed distance function (SDF)
sdf = sRectangle(0,W,0,H);
 
%% generate mesh
msh = Mesh(sdf,'Quads',1);
msh = msh.generate();

%% converting mesh to fem
fem = Fem(msh,'TimeStep',1/20,'SolverPlot',false,'ShowProcess',0);
      
%% assign boundary conditions
fem.AddConstraint('Support',fem.FindNodes('Bottom'),[0,1]);
fem.AddConstraint('Support',fem.FindNodes('Left'),[1,0]);
fem.AddConstraint('Displace',fem.FindNodes('Top'),[0,dL]);

%% add output
fem.AddConstraint('Output',fem.FindNodes('NW'),[0,0]);

%% material fitting
p0 = [0.1,0.01];

low = [1e-3,1e-3];
upp = [inf,inf];

P   = fmincon(@(p) Objective(fem,p,Data.Svm),p0,...
    [],[],[],[],low,upp);

% P   = fminunc(@(p) Objective(fem,p,Data.Svm),p0);

fprintf('Found Yeoh parameters \n');
fprintf('C1 = %3f3\n',P(1));
fprintf('C2 = %3f3\n',P(2));

%% plot optimization results
[lam,svm] = solverJob(fem,P);

figure(101); clf;
subplot(1,2,1);
fem.show('Uy'); 
subplot(1,2,2);

plot(lam,svm,'-','linewidth',2,'Color',col(1)); hold on;
plot(lam,Data.Svm,'.','MarkerSize',15,'Color',rgb2gray(col(4))); 

xaxis('Elongation strain $\lambda$','-')
yaxis('Engineering stress','Mpa');
axis([1,6,0,25]);

legend('Uni-axial Exp.','Material fit','Fontsize',16,...
    'Interpreter','Latex','Location','NorthWest');
grid on;

%% cost function
function J = Objective(fem,p,SvmExp)
% solve fem with parameters p
[~,svm] = solverJob(fem,p);

% compute cost function LEAST-SQUARES -> min J:= \SUM ||Sexp - S||_2
E = mean(SvmExp(:,1) - svm,1).';
J = E.'*E;
% print cost
fprintf('Cost function = %2f3 \n',J);
end

%% solver job routing
function [lam, svm] = solverJob(fem,p)
% assign material
fem.Material = YeohMaterial('C1',p(1),'C2',0.01);
% solving
fem.solve(); 
% extract
lam = (fem.BdBox(4)+fem.Log.Uy)/fem.BdBox(4);
svm = fem.Log.Svm;
end
