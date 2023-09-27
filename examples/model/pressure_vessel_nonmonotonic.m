clr;

fld = Fluidics();
fld.system.Control = @Controller;
fld.solver.TimeHorizon = 10;
fld.solver.TimeStep = 1/30;

R = 30; % mm
V = @(x) (4/3) * pi * R^3 * (1 + x^3 * tanh(x));
fld = fld.setPV(V);
fld = fld.setRegulator('on','pressure');

while fld.solver.Time < fld.solver.TimeHorizon
    fld = updateStatesFlow(fld);
end

t = fld.solver.sol.tout;
y = fld.solver.sol.yout;

subplot(1,2,1); ezplot(V); 
subplot(1,2,2); plot(t,y(:,1));

function u = Controller(x)
    u(1) = 10;
end