clr;
%% SDF
% P = polystar(3);
% sdf = sPolyline(P);

sdf = sCircle(1);
%%

sdf.show();
V = sdf.Node;

[T,N,B] = sdf.normal(V);

% P = [0,1;-1,0];
% N = (P*T.').';
% 
quiver(V(:,1),V(:,2),T(:,1),T(:,2)); hold on;
%quiver(V(:,1),V(:,2),N(:,1),N(:,2));
quiver(V(:,1),V(:,2),N(:,1),N(:,2));