clr;
%%
T  = 60;
dt = 1/200;
x0 = [0.5 0.5];

%%
f = @(x,u,t) [x(2); -10*sin(x(1)) - cos(x(2))];
sys = StateSpace(f,x0);

%%
mdl = Model(sys,'TimeEnd',T,'TimeStep',dt);
mdl = mdl.simulate();

%%
[ts,ys] = ode23t(@ODE,0:dt:T,x0);

%%
figure(103);
t = mdl.Log.t;
x = mdl.Log.x(:,1);

plot(t,x,'LineW',3); hold on;
plot(ts,ys(:,1),'--','LineW',3); hold on;

function dx = ODE(t,x)
 dx(1,1) = x(2);
 dx(2,1) = -10*sin(x(1)) - cos(x(2));
end