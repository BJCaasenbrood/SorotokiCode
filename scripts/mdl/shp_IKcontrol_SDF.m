clr; beep off;
%% settings
L = 1;   % manipulator length
M = 20;  % number of modes
N = 101; % grid on SR

%% build basis
x = linspace(0,1,N).';
Y = GenerateFunctionSpace(x,N,M,L);

%% desired SE3
Sd = sCircle(0.1,-0.3,0.1);  % desired enveloping SDF
sp = sSphere(0.1,0,0.3,0.1); % offset due to occupance of soft arm

%% soft sobotics shapes
figure(101); subplot(1,2,1);
shp = Shapes(Y,[0,M,0,0,0,0]);      % generate basis
rig = setupRig(M,L,[0,M,0,0,0,0]);  % rig the soft arm

q = 1e-4*sort(rand(shp.NDim,1));
%q(2) = -2;
e = 0.1;

%% solve IK
obj = Gmodel(sp);
obj.Texture = diffuse(0.925);
obj.bake.render();

BdBox = [-0.1,0.4,-0.2,0.2,0,0.6];

h = [];
EE = [];
k  = 0;

while norm(e) > 1e-3 && k < 400
    
    % update iteration
    k = k + 1;
    
    % compute Cosserat configuration
    [g, J] = shp.string(q);
    
    % extract positions
    p = reshape(g(1:3,4,:),3,[]).';
    
    % compute closest points
    V = Sd.Node;
    [XY,D]  = ClosestPointOnSDF(Sd,p(:,[1 3]));
    [T,N,B] = Sd.normal(XY);
    
    % update IK control law
    E(k) = 0;
    dq   = 0;
    
    A = [];
    b = [];
    
    for ii = 1:shp.NNode
        
        pd = [XY(ii,1);0;XY(ii,2)];
        Rd = [T(ii,1),B(ii,1),N(ii,1);
              T(ii,3),B(ii,3),N(ii,3);
              T(ii,2),B(ii,2),N(ii,2)];
        
        gd = SE3((Rd),pd);
        
        [dr, dE] = EnergyController(...
            g(:,:,ii), gd,...
            J(:,:,ii));                     
        
        dq = dq + dr;
        
        E(k) = E(k) + (dE.'*dE);
    end
    
    [Et,R] = shp.tangentPoint(q,...
        [XY(:,1),0*XY(:,1),XY(:,2)]);
    
    subplot(1,2,1);
    rig = rig.computeFK(q);
    rig = rig.update();
    
    % setup figure
    setupFigure(BdBox);
    
    hold on;
    if isempty(h)
        h = plot3(p(:,1),p(:,2),-p(:,3),'b-','LineW',3);
    else
        set(h,'XData',p(:,1));
        set(h,'YData',p(:,2));
        set(h,'ZData',-p(:,3));
    end
    
    subplot(1,2,2);
    plot(shp.Sigma,R,'LineW',3);
    axis([0 1 0 0.25]);
    drawnow;
    
    % compute update state and compute error
    q = q + dq;
    e = norm(abs(dq));
     
end

function [dq, E] = EnergyController(g,gd,J)
    
    k1   = 0.002;
    k2   = 0.075;
    lam1 = 1;
    
    % conditioner
    W  = 1;
    Kp = diag([k1,k1,k1,k2,k2,k2]);

    Xi = logmapSE3(g\gd);
    Fu = Kp*tmapSE3(W*Xi)*isomse3(W*Xi);

    dq = lam1*J.'*Fu;
    
    E = Kp*isomse3(Xi);
    
end

function [XY, D] = ClosestPointOnSDF(sdf,P)
    [XY, D] = distance2curve(sdf.Node,P);
end

function setupFigure(B)
    axis equal;
    axis(B);
    background(gitpage);
end

function rig = setupRig(M,L,Modes)

gmdl = Gmodel('Arm.stl');
gmdl.Alpha = 0.5;
 
N = 100;
X = linspace(0,L,N)';
Y = GenerateFunctionSpace(X,N,M,L);

shp = Shapes(Y,Modes,'L0',L);
 
rig = Rig(@(x) shp.string(x),'Domain',L);

rig = rig.add(gmdl);
rig = rig.parent(1,0,0);
rig = rig.parent(1,1,1);

rig    = rig.texture(1,1.1*base);
rig.g0 = SE3(roty(pi/2),zeros(3,1));

rig = rig.render();

view(20,15);

end

function Y = GenerateFunctionSpace(X,N,M,L)
% loop over functional space
Y = zeros(N,M);

for ii = 1:M
   %Y(:,ii) = chebyshev(X/L,ii-1); % chebyshev
   Y(:,ii) = pcc(X/L,ii,M); %chebyshev(X/L,ii-1); % chebyshev
end

% ensure its orthonormal (gram–schmidt)
Y = gsogpoly(Y,X);
end

