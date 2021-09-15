%
clr;
%% load object
obj = Gmodel('dzhanibekov.stl','Shading','Vertex');
obj.Texture = metal;

obj.bake(); obj.render();

%% bbl
[t,y] = ode45(@ODE,linspace(0,50,5e3),[0 1e-5 1]);  

figure(103);
plot(t,y);

%% plotting
R  = eye(3);
dt = mean(diff(t));

for ii = 1:fps(t,15):length(t)
    
    % update law rotation matrix
    R = R*expm(dt*so3(y(ii,:)));
    
    obj.reset();
    obj = Blender(obj,'SO3',R);
    
    obj.update();
    drawnow();
    
end


function dx = ODE(t,x)
J = 1e-6*diag([0.2,0.01,0.1]);

dx(1:3,1) = -J\(so3(x)*J*x);
end