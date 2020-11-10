clr;
%% preview
obj = Gmodel('Arm.stl');

%% set texture
obj.Texture = grey;
obj.bake().render().update();

%% transform
X = linspace(0,1,1e3);
gse3 = spiral(2,X(:));
id = knnsearch(X(:),obj.Node(:,3));
obj = Blender(obj,'Sweep', {id,gse3});

obj.update();
axis tight;

function g = spiral(a,t)
t = a*t;
g = zeros(length(t),7);
g(:,5:end) = [0.5*t,0.05*sqrt(t).*sin(pi*t),0.05*sqrt(t).*cos(pi*t)];
end