clr
%%
obj = Gmodel('Bunny.stl');
obj.Alpha = 0.5;

obj.bake.render();

F = obj.Element;
V = obj.Node;

B = obj.get('BdBox');

N = 4e3;
Pts = [];
while size(Pts,1) < N
    P = [random(B(1),B(2),10).',...
         random(B(3),B(4),10).',...
         random(B(5),B(6),10).'];
    h = inpolyhedron(F,V,P);

    Pts = [Pts; P(h,:)];

    size(Pts,1)
end

hold on;
fplot(Pts,'bo');

%%
clf;
SDF = @(x) dMetasphere(x,Pts(:,1),Pts(:,2),Pts(:,3),10 + Pts(:,3)*0,1);

sdf = Sdf(SDF,'BdBox',B,'Quality',60);
sdf.show();