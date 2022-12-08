clr; warning off;
%% settings
M = 7; 
N = 70;   % grid on SR

%% build basis
x = linspace(0,1,N).';
Y = zeros(N,M);

for ii = 1:M
    Y(:,ii) = x.^(ii-1);          % affine curvature
end

%% desired SE3
Sd = sCircle(0.25,-0.1,0.075);

gd = SE3(roty(pi/2),[0.5,0,-0.2]);
pd = gd(1:3,4);
%% soft sobotics shapes
shp = Shapes(Y,[0,M,0,0,0,0]);
q   = 0.01*ones(shp.NJoint,1);
e   = 0.1;

%% solve IK
figure(103); BdBox = [0,1.2,-.75,.75];

while norm(e) > 1e-3
    cla;
    
    % compute Cosserat configuration
    [g,J] = shp.string(q);
    
    % extract positions
    p  = reshape(g(1:3,4,:),3,[]).';
    
    % compute SDF tangent mapping
    [X,Y,Z] = TangentMap(p,BdBox);
    
    % compute closest points
    [XY,D]  = ClosestPointOnSDF(Sd,p(:,[1 3]));
    [T,N,B] = Sd.normal(XY);
    
    % plotting
    cplane(X,Y,Z); hold on;
    plot(p(:,1),p(:,3),'w-','LineW',4); 
    plot(pd(1),pd(3),'wo','LineW',2,'MarkerS',10);
   
    % setup figure
    setupFigure(BdBox);
    
    % update IK control law
    [dq, E] = EnergyController(g(:,:,end),gd,J(:,:,end));
    q = q + dq;
    
    % compute error
    e = abs(0.1*e - norm(abs(E)));

end

function [dq,E] = EnergyController(g,gd,J)
    k1 = 1;
    k2 = 0.1;
    lam1 = 1;
    lam2 = 1e-3;
    Kp = diag([k1,k1,k1,k2,k2,k2]);

    Xi = logmapSE3(g\gd);
    Fu = Kp*tmapSE3(Xi)*wedge(Xi);

    dq = lam1*J.'*((J*J.' + lam2*eye(6))\Fu);
    
    E = Kp*wedge(Xi);
end

function [X,Y,T] = TangentMap(p,BdBox)
    n = 150;
    x = linspace(BdBox(1),BdBox(2),n);
    y = linspace(BdBox(3),BdBox(4),n);

    [X,Y] = meshgrid(x,y);
    Pv = [X(:),Y(:)];

    Sf = sPolyline(p(:,[1 3]));
    Sf.BdBox = BdBox;

    [~,~,~,Z] = Sf.normal(Pv);
    T = reshape(Z,n,n);
end

function [XY,D] = ClosestPointOnSDF(sdf,P)
[XY,D,~] = distance2curve(sdf.Node,P);
end

function setupFigure(B)
    axis equal;
    axis(B);
    drawnow;
    colormap(turbo(0));
    pause(0.1);
end
