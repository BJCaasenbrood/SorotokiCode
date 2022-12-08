clr;
% CODE: Example 2 -- Task-space controller
% assign desired setpoint
Xd  = [30;10;10]; 

% build continuum shape
POD = chebyspace(60, 3);            % POD basis
shp = Shapes(POD,[0,3,3,0,0,0], ... % pure bending XY
                'Length', 100, ...
                'Material', Dragonskin10);

% geometry and boundary conditions
shp = shp.setRadius(5);   % R = 8 mm
shp = shp.setRamp(0.75);  % R is 3/4 at s=L
shp = shp.addGravity();

% model composer
mdl = Model(shp,'TimeStep', 1 / 60, ...
                'TimeEnd',  5, ...
                'Controller', @(x) tau(x, Xd) );

% magic ;)
mdl = mdl.simulate();

for ii = 1:numel(mdl.Log.t)  
    shp = shp.render(mdl.Log.x(ii,1:6));
    
    if ii == 1
        plotpoint(Xd); 
    end
    
    drawnow;
    view(30,30);
    axis([0 100 -5, 5, -10, 30]);
    
    if ii == 1
        gif('test.gif','frame',gcf,'nodither','DelayTime',1/60);
    else
        gif;
    end
    
    
end

% task-space controller (called by solver)
function tau = tau(mdl,Xd)
  log = mdl.Systems{1}.Log;
  X   = log.FK.g(1:3,4,end);
  Jv  = log.FK.J(4:6,:,end);
  Vq  = log.PH.dVdq;

  tau = Jv.'*(1e-3 * (Xd(:) - X) - 1e-4 * Jv * log.dq) + Vq;
end
