clr; warning off;
%% settings
m  = 3;    % number of modes
n  = 250;  % image grid
N  = 30;   % grid on SR
R  = 0.1;  % radius
dr = 2/n;  % delta r
x0 = 0;    % x0 of circle
y0 = -0.4;  % y0 of circle

%% build basis
x = linspace(0,1,N).';
Y = [];

for ii = 1:3
    Y(:,ii) = legendre(x,ii-1);   % legendre 
end

%% Soft Robotics Shapes
shp = Shapes(Y,[0,3,0,0,0,0]);

%% SDF of object
sdf = sCircle(x0,y0,R);
S = sdf.Node;
x = linspace(0,1,size(S,1));

%% evaluate geometric distance
[P,X,Y,I,D,Z,U,T] = GeometricSR(shp,sdf,[3,2,0],S);

% error
E = (Z-U).^2;

% interpolate
%I  = griddedInterpolant(X',Y',reshape(I,n,n)');
%I  = griddedInterpolant(P,I);
D  = griddedInterpolant(X',Y',reshape(D,n,n)');
Tx = griddedInterpolant(X',Y',reshape(T(:,1),n,n)');
Ty = griddedInterpolant(X',Y',reshape(T(:,2),n,n)');
Z  = griddedInterpolant(X',Y',reshape(Z,n,n)');
U  = griddedInterpolant(X',Y',reshape(U,n,n)');
E  = griddedInterpolant(X',Y',reshape(E,n,n)');

Sp = interp1(x,S,I);

%% plotting
% figure(101); clf;
% subplot(1,3,1);
% cplane(X,Y,D(X,Y)); hold on; 
% PlotSR(P,S); axis equal; axis tight;
% plot(Sp(:,1),Sp(:,2),'r-');
% 
% subplot(1,3,2);
% cplane(X,Y,U(X,Y)); hold on; 
% PlotSR(P,S); axis equal; axis tight;



%% function
function PlotSR(P,S)
    plot(P(:,1),P(:,2),'-w','LineW',5); hold on;
    plot(P(end,1),P(end,2),'w.','MarkerSize',30); 
    plot(P(1,1),P(1,2),'w.','MarkerSize',30); 
    
    plot(S(:,1),S(:,2),'w','LineW',3); 
    axis equal tight
    colormap(turbo);
end

%% 
function [P,X,Y,I,D,Z,U,T] = GeometricSR(shp,sdf,Q,S)
g = shp.string(Q); 
p = reshape(g(1:3,4,:),3,[]).';
P = [p(:,1),p(:,3)];

BdBox = [-0.25 1.25 -1.5 1.5];

Sf = sPolyline(P);
sdf.BdBox = BdBox; 
Sf.BdBox  = BdBox;

% reconstruct theta-field
n = 250;
x = linspace(BdBox(1),BdBox(2),n);
y = linspace(BdBox(3),BdBox(4),n);

[X,Y] = meshgrid(x,y);
Pv = [X(:),Y(:)];

[T,Z] = Sf.tangent(Pv);
[~,U] = sdf.tangent(Pv);

D = Sf.eval(Pv); 
D = D(:,end);

[~,~,I] = distance2curve(S,P);
end

