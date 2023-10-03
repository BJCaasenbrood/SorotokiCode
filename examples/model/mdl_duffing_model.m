clr;
% parameters
[amp, b, alpha, beta, w] = deal(0.45, 0.5, -1, 1, 2);

% flow and hessian function
f = @(x,u,t) [ x(2); -b*x(2)-alpha*x(1)-beta*x(1)^3+amp*sin(w*t)];
H = @(x,u,t) [ 0, 1; -alpha-3*beta*x(1)^2, -b];

% build nonlinear statespace model \lambda
sys = StateSpace({f,H},[0.5021; 0.17606]);
mdl = Model(sys);

% setting variables
mdl = mdl.set('MaxIteration',2,...
        'RelTolerance',1e-9,...
        'TimeStep',0.1,...
        'TimeHorizon',50);

mdl = mdl.simulate();

% plotting duffing solution
X = mdl.solver.sol.yout(:,1);
Y = mdl.solver.sol.yout(:,2);
plot(X,Y);
