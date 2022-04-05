clr; beep off;
%% settings
L = 1;   % manipulator length
M = 4;   % number of modes
N = M*8; % grid on SR

%% build basis
x = linspace(0,1,N).';
Y = GenerateFunctionSpace(x,N,M,L);

%% desired SE3
Sd = sCircle(0.1,-0.3,0.15);  % desired enveloping SDF
sp = sSphere(0.1,0,0.3,0.09); % offset due to occupance of soft arm

%% soft sobotics shapes
shp = Shapes(Y,[0,M,0,0,0,0]);      % generate basis
rig = setupRig(M,L,[0,M,0,0,0,0]);  % rig the soft arm

q = 1e-4*sort(rand(shp.NDim,1));
e = 0.1;

%% solve IK
figure(101); 

obj = Gmodel(sp);
obj.Texture = diffuse(0.925);
obj.bake.render();

BdBox = [0.3,0.6,-0.4,0.4,0,0.6];

EE = [];
k  = 0;

while norm(e) > 1e-3 && k < 400
    %clf;
    
    % update iteration
    k = k + 1;
    
    % compute Cosserat configuration
    [g, J] = shp.string(q);
    
    % extract positions
    p  = reshape(g(1:3,4,:),3,[]).';
    
    % compute closest points
    V = Sd.Node;
    [XY,D]  = ClosestPointOnSDF(Sd,p(:,[1 3]));
    [T,N,B] = Sd.normal(XY);
    
    % plotting
    %plotSE2(g);
    %Sd.showcontour();
    
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
            g(:,:,ii),...
            gd,...
            J(:,:,ii));                     
        
        dq = dq + dr;
        
        E(k) = E(k) + (dE.'*dE);
    end
   
    q = q + dq;
    
    rig = rig.computeFK(q);
    rig = rig.update();
    
    % compute update state and compute error
    q = q + dq;
    e = norm(abs(dq));
    
    % setup figure
    setupFigure(BdBox);
    drawnow;
    
    if k == 1, gif('shp_ik_sdf.gif','frame',gcf,'nodither');
    else, gif;
    end
     
end

function [dq, E] = EnergyController(g,gd,J)
    
    k1   = 0.002;
    k2   = 0.1;
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
 
N = 200;
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
   Y(:,ii) = chebyshev(X/L,ii-1); % chebyshev
end

% ensure its orthonormal (gram–schmidt)
Y = gsogpoly(Y,X);
end

