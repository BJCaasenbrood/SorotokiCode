clr;
%% settings
N = 120;
R = 0.15;
x0 = 0; 
y0 = 0.5;

%% build basis
x = linspace(0,1,N).';
Y = [];

for ii = 1:4
    %Y(:,ii) = cos(2*pi*(ii-1)*x); % sinusoidal
    %Y(:,ii) = legendre(x,ii-1);   % legendre
    Y(:,ii) = chebyshev(x,ii-1);  % chebyshev
    %Y(:,ii) = x.^(ii-1);           % affine curvature
end

Y = gsogpoly(Y);

%% Shapes
Table = [0,4,0,0,0,0];
shp = Shapes(Y,Table);

%% 
g  = shp.string([2,-3,0,2]); 
th = 2*pi*x - pi/2;

p = reshape(g(1:3,4,:),3,[]).';

P = [p(:,1),p(:,3)];
S = [R*cos(th)+x0,R*sin(th)+y0];
%% 
BDX = [-0.4 .4]*2-.4;
BDY = -0.75 + [0 1]*(BDX(2)-BDX(1))*1.41/2;
BdBox = [BDX BDY];

Cf = sCircle(x0,y0,R);
Sf = sPolyline(P);
Cf.BdBox = BdBox; 
Sf.BdBox = BdBox;

%% reconstruct theta-field
n = 500;
x = linspace(BdBox(1),BdBox(2),n);
y = linspace(BdBox(3),BdBox(4),n);

[X,Y] = meshgrid(x,y);
Pv = [X(:),Y(:)];

tic
[~,Z] = Sf.tangent(Pv);
[~,U] = Cf.tangent(Pv);

SDF = Cf;

D = Sf.eval(Pv); D = D(:,end);
%U = Cf.eval(Pv); U = U(:,end);
toc

[~,~,I] = distance2curve(P,Pv);

%% plotting
figure(101); clf;
% subplot(1,4,1);
% %E = (U-Z);
% % E(Z>0) = 0;
% %E = U;
% cplane(X,Y,reshape(I,n,n)); hold on; 
% PlotSR(P,S);
% %caxis([-pi,pi]); 
% %caxis([min(U),0]); 
% caxis([0,1]);
% colormap(turbo());
% axis equal; axis tight; axis off 
% colorbar('location','southoutside');
% title('\sigma-field','fontsize',18);
% 
% subplot(1,4,2);
% E = I.*(U-Z).^2;
% E(D<0) = 0;
% cplane((X),(Y),reshape(E,n,n)); hold on; 
% PlotSR(P,S);
% caxis([0,max(E)]); 
% %caxis([min(Z),0]); 
% colormap(turbo());
% axis equal; axis tight; axis off 
% colorbar('location','southoutside');
% title('\Theta-field','fontsize',18);
% 
% subplot(1,4,3);
% cplane((X),(Y),reshape((U),n,n)); hold on; 
% PlotSR(P,S);
% caxis([-pi,pi]); 
% %caxis([min(Z),0]); 
% colormap(turbo());
% axis equal; axis tight; axis off 
% colorbar('location','southoutside');
% title('SDF-field','fontsize',18);
% 
% 
% subplot(1,4,4);
cplane((X),(Y),reshape((Z),n,n)); hold on; 
PlotSR(P,S);
caxis([-pi,pi]); 
%caxis([min(Z),0]); 
colormap(turbo());
axis equal; axis tight; axis off 
%colorbar('location','southoutside');
title('SDF-field','fontsize',18);


%% function
function PlotSR(P,S)
    plot(P(:,1),P(:,2),'-w','LineW',5); hold on;
    plot(P(end,1),P(end,2),'w.','MarkerSize',30); 
    plot(P(1,1),P(1,2),'w.','MarkerSize',30); 
    
    %plot(S(:,1),S(:,2),'w','LineW',3); 
end

%%
% N = repmat([0,1],[length(Pv),1]);
% %Z = acos(dot(T',N').');

%surf(X,Y,reshape(Z,n,n),'Linestyle','none');
%C1 = sCircle(0.1,-0.4,0.2);
%C1.showcontour();

%%
% th = RotationAngle(R);
% 
% N  = 250;
% Pp = [p(:,1),p(:,3)];
% 
% X = linspace(-0.5,1.5,N);
% Y = linspace(-.75,0.75,N);
% [XX,YY] = meshgrid(X,Y);
% 
% Pv = [XX(:),YY(:)];
% 
% [xy,distance,t_a] = distance2curve(Pp, Pv);
% 
% id = round(t_a*(shp.NNode-1))+1;
% 
% relDist = Pv - Pp(id,:);
% relTh   = th(id,:);
% 
% Z = RelativeRotation(relTh,relDist);
% ZZ = reshape(Z,[N,N]);
% caxis([-pi,pi])
% surf(XX,YY,-ZZ,'linestyle','none')
% 
% colormap(turbo)
% view(0,90);
% 
% function theta = RotationAngle(Rot)
% N  = length(Rot);
% 
% for ii = 1:N
%     T = logm(Rot{ii});
%     theta(ii,1) = T(3,1);
% end
% 
% end
% 
% function Rp = RelativeRotation(theta,rd)
% Rp = theta*0;
% for ii = 1:length(theta)
%     th = theta(ii);
%     RT = [cos(th),sin(th);-sin(th),cos(th)];
%     D = RT.'*rd(ii,:).';
%     Rp(ii,1) = atan2(D(1),D(2));
% end
% 
% %Rp = wrap(Rp,[-pi,pi]);
% 
% end
% 
% 
