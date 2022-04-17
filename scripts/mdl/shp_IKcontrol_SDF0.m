clr; beep off; 
%%  settings
M = 7;   % number of modes
N = 50;  % grid on SR

%% build basis
x = linspace(0,1,N).';
Y = zeros(N,M);

for ii = 1:M
    Y(:,ii) = chebyshev(x,ii-1); % legendre 
end

%% desired SE3
sdf = sCircle(0.25,-0.1,0.075);

%% soft sobotics shapes
shp = Shapes(Y,[0,M,0,0,0,0]);

lam = 0;
q  = 1e-4*sort(rand(shp.NDim,1));
e  = 0.1;
%% 
zeta = 0.02;
Kee  = shp.Ktt;
Cmat = zeros(shp.NDim,shp.NDim);

for ii = 1:shp.NNode
   N = shp.Ba*shp.Theta(shp.Sigma(ii));
   Cmat = Cmat + double(zeta*N.'*Kee*N*(1/shp.NNode));
end

%% kinmatic solver
figure(101);

for ii = 1:100
    
    % eval shape
    [g,J] = shp.string(q(:));
    
    % get projection on SDF
    [XY] = ClosestProjection(sdf,g);
    
    % get normal field
    [T,N,B] = sdf.normal(XY);
    
    % construct desired manifold
    gd = ConstructDesiredSE3(XY,T,N,B);
    
    c = 1;
    k = 1;
    
    while abs(c) > 1e-3
        
        % eval flow
        [c,dcdq] = Constraint(q,sdf,shp);
        
        % eval flow
        [f,dfdq] = Flow(g,gd,J);
        
        A = [Cmat,-dcdq.';-dcdq,1e-12];
        b = [f*smoothstep(0.02*ii)-0*dcdq.'*lam;c];
        
        %clc[L,D,P] = ldl(A,'vector');
        %     %minDiag = full(min(diag(D)));
        %dx = sparse(P,1,(L'\(D\(L\b(P)))));
        
        dx = A\b;
        
        lam = 0.05*dx(end);
        q   = q + dx(1:shp.NDim);
        k   = k + 1;
    end

    cla;
    plotSE2(g); hold on;
    sdf.showcontour();
    
    axis equal;
    axis([-1,1,-1,1]);
    drawnow;
   
end


%% functions
function [c] = Cfunc(q,sdf,shp)
    % stiffness
    k = 0.0;
    
    % compute string
    [g] = shp.string(q(:));
    
    % extract positions
    p = reshape(g(1:3,4,:),3,[]).';
    N = reshape(g(1:3,3,:),3,[]).';
    
    % eval distance
    d = sdf.eval([p(:,1),p(:,3)]);
    [T,~] = sdf.normal([p(:,1),p(:,3)]);
    
    d = d(:,end);
    
    % compute herzian contact force
    f = k*min(d,0);
    
    % compute equality constr.
    c = sum(dot(-(N.*f).',T.'));
end

function [c,dcdq] = Constraint(q,sdf,shp)
    
    de = 1e-3;
    
    c = Cfunc(q,sdf,shp);
    dq    = eps*eye(shp.NDim);
    dcdq  = zeros(shp.NDim,1);
    
    for ii = 1:shp.NDim
        ci = Cfunc(q(:) + dq(:,ii),sdf,shp);
        dcdq(ii) = (ci - c)/de;
    end
    
    dcdq = dcdq.';

end

function [XY] = ClosestProjection(sdf,g)

% extract positions
p = reshape(g(1:3,4,:),3,[]).';

% compute closest points
[XY,~] = distance2curve(sdf.Node,p(:,[1 3]));

end

function gd = ConstructDesiredSE3(XY,T,N,B)

gd = zeros(4,4,size(XY,1));

for ii = 1:size(XY,1)
    
        % get desired distance
        pd = [XY(ii,1);0;XY(ii,2)];
    
        Rd = [T(ii,1),B(ii,1),N(ii,1);
              T(ii,3),B(ii,3),N(ii,3);
              T(ii,2),B(ii,2),N(ii,2)];
    
        % get desired SE(3)
        gd(:,:,ii) = SE3(Rd,pd);
   
end
    
end

function [f,dfdq] = Flow(g,gd,J)

lam1 = 0.25;
lam2 = .01;

f    = zeros(size(J,2),1);
dfdq = zeros(size(J,2),size(J,2));

for ii = 1:size(g,3)
    
    % compute force residual
    dr = EnergyBasedController(g(:,:,ii),gd(:,:,ii));
    
    df = lam1*J(:,:,ii); %
    %df = lam1*J(:,:,ii).'*inv((J(:,:,ii)*J(:,:,ii).' + lam2*eye(6)));
    
    %f    = [f;dr];
    %dfdq = [dfdq;df];
    f    = f + J(:,:,ii).'*dr;
    dfdq = dfdq + J(:,:,ii).'*J(:,:,ii);
    
    % update dq
    %dq = dq + dr;
end

end

function [dr] = EnergyBasedController(g,gd)
    
    k1 = 0.01;
    k2 = 1.5;
    
    Kp = diag([k1,k1,k1,k2,k2,k2]);

    Xi = logmapSE3(g\gd);
    dr = Kp*tmapSE3(Xi)*isomse3(Xi);

end