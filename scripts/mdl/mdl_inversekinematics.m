clr; beep off;
%% settings
M = 4;    % number of modes
N = M*8;   % grid on SR

%% build basis
x = linspace(0,1,N).';
Y = zeros(N,M);

for ii = 1:M
    Y(:,ii) = chebyshev(x,ii-1); % legendre 
    %Y(:,ii) = pcc(x,ii,M);       % piece-wise constant
end

Y = gsogpoly(Y);

%% desired SE3
Sd = sCircle(0.25,-0.1,0.075);

%% soft sobotics shapes
shp = Shapes(Y,[0,M,0,0,0,0]);

q  = 1e-4*sort(rand(shp.NDim,1));
e  = 0.1;

%% solve IK
figure(103); BdBox = [0,1.2,-.75,.75];

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
    subplot(1,2,1); hold on;
    %cplane(X,Y,Z);  hold on;
    %plot(p(:,1),p(:,3),'k-','LineW',2);
    %plot(p(Pts,1),p(Pts,3),'k.','MarkerS',15);
    plotSE2(g);
    plot(V(:,1),V(:,2),'k--','MarkerS',10);
    
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
    
    subplot(1,2,2);
    plot(E,'LineW',3);
    %axis equal;
    axis([0 100 0 0.05]);
    title('Energy difference');
    drawnow;
   
end

function [dq, E] = EnergyController(g,gd,J,k)
    
    k1 = 0.000;
    k2 = 0.25;
    lam1 = 1;
    
    % conditioner
    W  = smoothstep(k/5+0.1);
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
end
