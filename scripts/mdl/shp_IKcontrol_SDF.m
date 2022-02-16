clr; beep off;
%% settings
M = 5;    % number of modes
N = 60;   % grid on SR

%% build basis
x = linspace(0,1,N).';
Y = zeros(N,M);

for ii = 1:M
    Y(:,ii) = chebyshev(x,ii-1); % legendre 
end

Y = gsogpoly(Y);

%% desired SE3
Sd = sCircle(0.5,-0.2,0.075);

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
    [g,J] = shp.string(q);
    
    % extract positions
    p  = reshape(g(1:3,4,:),3,[]).';
    
    % compute SDF tangent mapping
    %[X,Y,Z,Sp] = TangentMap(p,BdBox);
    
    % compute closest points
    V = Sd.Node;
    [XY,D]  = ClosestPointOnSDF(Sd,p(:,[1 3]));
    [T,N,B] = Sd.normal(XY);
    
    % plotting
    subplot(1,2,1); hold on;
    %cplane(X,Y,Z);  hold on;
    plot(p(:,1),p(:,3),'k-','LineW',2);
    plot(V(:,1),V(:,2),'k--','MarkerS',10);
    colormap(bluesea(0));
    
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
        
        %dist = Sd.eval([g(1,4,ii),g(3,4,ii)]); 
        %dist = clamp(-dist(end),0,Inf);  
        %Fd = [0;0;0;Rd*[0;0;dist]];
        
        %Null = (eye(7) - pinv(J(:,:,ii))*J(:,:,ii));
                                
        %lambda = (J(:,:,ii).'*Fd);
        dq = dq + dr ;%+ lambda;
        
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
    %colorbar('location','NorthOutside');
    
    subplot(1,2,2);
    plot(E,'LineW',3);
    axis equal;
    axis([0 100 0 45]);
    title('Energy difference');
    drawnow;
    
%     if k == 1
%         gif('grasping.gif','frame',gcf,'nodither');
%         background('w');
%     else
%        gif; 
%     end
    
end

function [dq,E] = EnergyController(g,gd,J,k)
    
    k1 = 0.01;
    k2 = 1;
    
    lam1 = .3;
    lam2 = 1;
    
    Kp = diag([k1,k1,k1,k2,k2,k2]);

    Xi = logmapSE3(g\gd);
    Fu = Kp*tmapSE3(smoothstep(k/30+0.1)*Xi)*isomse3(smoothstep(k/30+0.1)*Xi);

    %dq = lam1*J.'*((J*J.' + lam2*eye(6))\Fu);
    
    dq = lam1*J.'*Fu;
    
    E = Kp*isomse3(Xi);
end

function [X,Y,T,Sf] = TangentMap(p,BdBox)
    n = 2;
    x = linspace(BdBox(1),BdBox(2),n);
    y = linspace(BdBox(3),BdBox(4),n);

    [X,Y] = meshgrid(x,y);
    Pv = [X(:),Y(:)];

    Sf = sPolyline(p(:,[1 3]));
    Sf.BdBox = BdBox;

    [~,~,~,Z] = Sf.normal(Pv);
    T = reshape(Z,n,n);
end

function [XY, D] = ClosestPointOnSDF(sdf,P)
    [XY, D] = distance2curve(sdf.Node,P);
end

function setupFigure(B)
    axis equal;
    axis(B);
    %drawnow;
    %colormap(bluesea(0));
    %pause(0.1);
end
