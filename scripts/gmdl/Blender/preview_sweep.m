clr;
%% preview
obj = Gmodel('SlenderRod.stl');

%% set texture
obj.Texture = prusa;
obj.bake().render().update();

%% transform
X = linspace(0,1,500);
gse3 = spiral(3,X(:));
id = knnsearch(X(:),obj.Node(:,3));
obj = Blender(obj,'Sweep', {id,gse3});

obj.update();
axis tight;

function g = spiral(a,t)
t = a*t;
g = zeros(length(t),7);
g(:,5:end) = [0.5*t,0.25*sqrt(t).*sin(pi*t),0.25*sqrt(t).*cos(pi*t)];
end