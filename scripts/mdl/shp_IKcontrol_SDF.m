clr; beep off;
%% settings
L = 1;   % manipulator length
M = 8;   % number of modes
N = 101; % grid on SR

%% build basis
x = linspace(0,1,N).';
Y = GenerateFunctionSpace(x,N,M,L);

%% desired SE3
R = 0.15;
D = 0.055;
Sd = sCircle(0.2,-0.15,R+D); % desired enveloping SDF
sp = sSphere(0.2,0,0.15,R);  % offset due to occupance of soft arm

%% soft sobotics shapes
figure(101); subplot(1,2,1);
shp = Shapes(Y,[0,M,0,0,0,0]);      % generate basis
rig = setupRig(M,L,[0,M,0,0,0,0]);  % rig the soft arm

q = 1e-4*sort(rand(shp.NDim,1));

%% solve IK
obj = Gmodel(sp,'Texture',diffuse(0.925));
obj.bake.render();

BdBox = [-0.1,0.4,-0.2,0.2,0,0.6];

h  = [];
EE = [];
k  = 0;
e  = 0.1;

while norm(e) > 1e-3 && k < 400
    
    % update iteration
    k = k + 1;
    
    % compute Cosserat configuration
    [g, J] = shp.string(q);
    
    % extract positions
    p = reshape(g(1:3,4,:),3,[]).';
    
    % compute closest points
    V = Sd.Node;
    
    [XY,~]  = ClosestPointOnSDF(Sd,p(:,[1 3]));
    [T,N,B] = Sd.normal(XY);
    
    % update IK control law
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
    end
    
    [Et,R] = shp.tangentPoint(q,...
        [XY(:,1),0*XY(:,1),XY(:,2)]);
    
    subplot(1,2,1);
    rig = rig.computeFK(q);
    rig = rig.update();
    
    V = rig.List{1}.Node;
    V = V(sp.intersect(V),:);
    
    [~,N] = sp.normal(V);
    
    delete(h); hold on;
    h = quiver3(V(:,1),V(:,2),V(:,3),...
                N(:,1),N(:,2),N(:,3),'b-');
    
    % setup figure
    setupFigure(BdBox);
    
    subplot(1,2,2);
    plot(shp.Sigma,R,'LineW',3);
    axis([0 1 0 0.25]);
    axis square; ax = gca;
    grid on;
    set(ax,'LineW',1.5);
    ax.YAxis.Exponent = -2;
    ax.XTickLabel = {'0',[],'L'};
    ax.FontSize = 12;
    drawnow;
    
    % compute update state and compute error
    q = q + dq;     
end

function [dq, E] = EnergyController(g,gd,J)
    
    k1   =0;
    k2   = 0.75;
    lam1 = 1;
    
    % conditioner
    W  = 1;
    Kp = diag([k1,k1,k1,k2,k2,k2]);

    Xi = logmapSE3(g\gd);
    Fu = Kp*tmapSE3(W*Xi)*wedge(W*Xi);

    dq = lam1*J.'*Fu;
    
    E = Kp*wedge(Xi);
    
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

gmdl = Gmodel('Arm.stl','ShowProcess',0,'Alpha',0.3);
 
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

