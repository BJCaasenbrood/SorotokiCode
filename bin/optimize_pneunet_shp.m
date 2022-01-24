clr;
%%
load('job.mat');

%%
x0 = [0,0,0];
q = [];
for ii = 1:251
N  = fem.Log.Node{ii};
R = fem.Log.Rotation{ii};

L = shp.get('Filter');
pd = L*N;
rd = {};

lst = 1:shp.Fem.NNode;
for jj = 1:shp.NNode
    
    RRe = 0;
    W   = L(jj,:);
    for kk = lst(abs(L(jj,:))>0)
        RRe = RRe + W(kk)*R{kk};
    end
    
    [Ur,~,Vr] = svd(RRe);
    Re = (Ur*Vr.');
    
    rd{jj,1} = Re;
end


figure(101); cla;
plot(N(:,1),N(:,2),'b.');
plot(pd(:,1),pd(:,2),'ro');

DiffStepSize = [1e-3,1e-5,1e-3];
%options = optimoptions('fminsearch','Display','iter');%,...%'Algorithm','quasi-newton',...
    %'FiniteDifferenceStepSize',DiffStepSize);

fun = @(x) Objective(x,shp,pd,rd,x0);
x   = fminsearch(fun,x0);

P = shp.string(x);
plot(P(:,1),P(:,3),'r-'); hold on;

x0 = x
axis equal;
pause(0.5);
q = [q;x(:).'];
end

function J = Objective(x,shp,pd,rd,x0)

[P, R] = shp.string(x);

W = 1 + 0*(shp.Sigma/max(shp.Sigma)).^2;

d1 = sqrt((P(:,1) - pd(:,1)).^2 + (P(:,3) - pd(:,2)).^2);
J1 = sum(W(:).*d1(:));

for ii = 1:shp.NNode
    Rdif = R{ii}*rd{ii}.';
    th(ii) = W(ii)*acos(0.5*trace(Rdif) - 0.5);
end

% d2 = sqrt((P(end,1) - rd(end,1)).^2 + (P(end,3) - rd(end,2)).^2);
% J1 = sum(d2);

J2 = sum(th);
Q = eye(numel(x0));
a = 0.5;
J = (1-a)*J1 + a*J2 + (x(:)-x0(:)).'*Q*(x(:)-x0(:));

%p = [P(:,1),P(:,3)];
%[~,dist] = distance2curve(p/1e3,rd/1e3,'linear'); 

%J = sum(abs(dist));

%d2 = sqrt((P(end,1) - rd(end,1)).^2 + (P(end,3) - rd(end,2)).^2);
%J = d2;
%IDX = knnsearch([P(:,1),P(:,3)],rd)

end