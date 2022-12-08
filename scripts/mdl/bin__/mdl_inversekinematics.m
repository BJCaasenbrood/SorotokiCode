clr; beep off;
%% settings
M = 4;    % number of modes
N = 10;   % grid on SR

%% build basis
x = linspace(0,1,N).';
Y = zeros(N,M);

for ii = 1:M
    Y(:,ii) = chebyshev(x,ii-1); % legendre 
end

Y = gsogpoly(Y);

%% desired SE3
Sd = sCircle(-0.25,-0.6,0.15);

%% soft sobotics shapes
Y = chebyspace(N,M);
shp = Shapes(Y,[0,M,0,0,0,0]);
shp.g0 = SE3(roty(pi/2),[0;0;0]);

q  = 1e-4*sort(rand(shp.NDim,1));
e  = 0.1;

%% solve IK
figure(103); BdBox = [-1,1,-1.5,.05];

EE = [];
k  = 0;

while norm(e) > 1e-3 && k < 100
    clf;
    
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
    %subplot(1,2,1); hold on;
    plotSE2(g,'xz',0.01);
    plot(V(:,1),V(:,2),'k','MarkerS',10);
    
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
        %plotSE2(gd);
        plot(pd(1),pd(3),'k.');
        
        [dr, dE] = EnergyController(g(:,:,ii),...
                                    gd,J(:,:,ii),k);                     
        
        dq = dq + dr;
        
        E(k) = E(k) + (dE.'*dE);
    end
   
    q = q + dq;
    
    % iteration counter
    fprintf(' iteration = %i \n', k);
    
    % compute update state and compute error
    q = q + dq;
    e = norm(abs(dq));
    
    % setup figure
    setupFigure(BdBox);
%     
%     subplot(1,2,2);
%     plot(E,'LineW',3);
%     
%     %axis equal;
%     axis([-50 50 0 0.05]);
%     title('Energy difference');
    drawnow;
   
end

function [dq, E] = EnergyController(g,gd,J,k)
    
    k1 = 0.00;
    k2 = 0.5;
    lam1 = 1e4;
    lam2 = 3e3;
    
    % conditioner
    W  = smoothstep(k/5+0.1);
    Kp = diag([k1,k1,k1,k2,k2,k2]);

    Xi = logmapSE3(g\gd);
    Fu = Kp*tmapSE3(W*Xi)*wedge(W*Xi);

    dq = lam1*J.'*((J*J.' + lam2*eye(6))\Fu);
    
    E = Kp*wedge(Xi);
    
end

function [XY, D] = ClosestPointOnSDF(sdf,P)
    [XY, D] = distance2curve(sdf.Node,P);
end

function setupFigure(B)
    axis equal;
    axis(B);
end