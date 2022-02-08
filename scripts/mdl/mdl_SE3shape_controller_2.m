clr; beep off;
%% settings
M = 10; 
N = 100;   % grid on SR

%% build basis
x = linspace(0,1,N).';
Y = zeros(N,M);

for ii = 1:M
    Y(:,ii) = chebyshev(x,ii-1); % legendre 
    %Y(:,ii) = x.^(ii-1); % legendre 
end

%Y = gsogpoly(Y);

%% desired SE3
Sd = sCircle(0.3,-0.1,0.075);

%% soft sobotics shapes
shp = Shapes(Y,[0,M,0,0,0,0]);
q   = 1e-6*rand(shp.NDim,1);
e   = 0.1;

%% optimization
figure(101);
x = fminunc(@(x) Objective(x,shp,Sd),q);

function Cost = Objective(q,shp,Sd)
    
     [g,J] = shp.string(q);
    
     p = reshape(g(1:3,4,:),3,[]).';
     
     [XY,D]  = ClosestPointOnSDF(Sd,p(:,[1 3]));
     [T,N,B] = Sd.normal(XY);
     
     F = 0;
     E = 0;

     for ii = 1:shp.NNode
         
         pd = [XY(ii,1);0;XY(ii,2)];
         Rd = [T(ii,1),B(ii,1),N(ii,1);
               T(ii,3),B(ii,3),N(ii,3);
               T(ii,2),B(ii,2),N(ii,2)];
         
         gd = SE3((Rd),pd);
         
         [dF, dE] = EnergyController(g(:,:,ii),...
             gd);
         
         %dq = dq + shp.Sigma(ii)*dr;
         F = F + (dF).'*dF;
         %E = E + dE.'*dE;
         %subplot(1,2,1); hold on;
         %plotSE2(gd,'xz',0.1);
     end
     
    % plotting
    %figure(101); 
    cla; hold on;
    %subplot(1,2,1); hold on;
    %cplane(X,Y,Z); hold on;
    %plot(p(:,1),p(:,3),'k-','LineW',4);
    plotSE2(g);
    plot(Sd.Node(:,1),Sd.Node(:,2),'k--','MarkerS',10);
    colormap(bluesea(0));
    axis equal; axis([-.15,1,-1,0.5]);
    drawnow;
    Cost = F
     
end

% %% solve IK
% figure(103); BdBox = [0,1.2,-.75,.75];
% 
% EE = [];
% k = 0;
% while norm(e) > 1e-3
%     clf;
%     
%     % update iteration
%     k = k + 1;
%     
%     % compute Cosserat configuration
%     [g,J] = shp.string(q);
%     
%     % extract positions
%     p  = reshape(g(1:3,4,:),3,[]).';
%     
%     % compute SDF tangent mapping
%     [X,Y,Z,Sp] = TangentMap(p,BdBox);
%     
%     % compute closest points
%     V = Sd.Node;
%     [XY,D]  = ClosestPointOnSDF(Sd,p(:,[1 3]));
%     [T,N,B] = Sd.normal(XY);
%     
%     % plotting
%     subplot(1,2,1); hold on;
%     cplane(X,Y,Z); hold on;
%     plot(p(:,1),p(:,3),'k-','LineW',4);
%     plot(V(:,1),V(:,2),'k--','MarkerS',10);
%     colormap(bluesea(0));
%     
%     % update IK control law
%     E(k) = 0;
%     dq = 0;
%     
%     A = [];
%     b = [];
%     
%     for ii = 1:shp.NNode
%         
%         pd = [XY(ii,1);0;XY(ii,2)];
%         Rd = [T(ii,1),B(ii,1),N(ii,1);
%               T(ii,3),B(ii,3),N(ii,3);
%               T(ii,2),B(ii,2),N(ii,2)];
%         
%         gd = SE3(flipud(Rd),pd);
%         
%         [dr, dE] = EnergyController(g(:,:,ii),...
%                                     gd,J(:,:,ii));                     
%         
%         dq = dq + shp.Sigma(ii)*dr;
%         
%         E(k) = E(k) + dE.'*dE;
%     end
%     
% %     dq = A.'*b;
% %     q_ = q;
% %     q = (1-0.01)*q_ + (0.01)*(q + dq);
%     
%     % iteration counter
%     fprintf(' iteration = %i \n', k);
%     
%     % compute update state and compute error
%     q = q + dq;
%     e = norm(abs(dq));
%     
%     % setup figure
%     setupFigure(BdBox);
%     
%     subplot(1,2,2);
%     plot(E,'LineW',3);
%     %axis equal;
%     axis([0 150 0 max(E)]);
%     drawnow;
%     
% end

function [Fu,E] = EnergyController(g,gd)
    
    k1 = 0.05;
    k2 = 1;
    
    Kp = diag([k1,k1,k1,k2,k2,k2]);

    Xi = logmapSE3(g\gd);
    Fu = Kp*tmapSE3(Xi)*isomse3(Xi);
    %Fu = isomse3(Xi);
    %dq = lam1*J.'*inv(J*J.' + lam2*eye(6))*Fu;
    %dq = lam1*J.'*Fu;
    
    E = Kp*isomse3(Xi);
end


function [XY, D] = ClosestPointOnSDF(sdf,P)
    [XY, D] = distance2curve(sdf.Node,P);
end

function setupFigure(B)
    axis equal;
    axis(B);
    %drawnow;
    colormap(bluesea(0));
    %pause(0.1);
end
