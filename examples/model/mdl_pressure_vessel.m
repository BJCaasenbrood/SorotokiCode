clr;
% set fluidics model
fld = Fluidics();
fld.system.Control = @Controller;
fld = fld.setRegulator('on','pressure');

fld.params.Kp = 0.1;

while fld.solver.Time < fld.solver.TimeHorizon
    fld = updateStatesFlow(fld);
end

t = fld.solver.sol.tout;
y = fld.solver.sol.yout;

fig;
plot(t,y(:,1));
ylim([0 15]);

function u = Controller(sys)
    t = sys.solver.Time;
    u = waveform(t,[0,10;5,-10;10,10]);
end