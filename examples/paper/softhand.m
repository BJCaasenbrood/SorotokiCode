clr; 
% load hand preset
[shp, obj] = preset.shapes.softhand;
f = figure(101);

dt = 1/30; % set initial timestep

% simulation loop
while true
    t = shp{1}.solver.Time;
    for ii = 1:5, 
      shp{ii}.system.fControl = [500; 0] * smoothstep(0.15 * t ); 
    end

    tic;
    shp = cellfun(@(x) x.update(dt), shp, 'UniformOutput', false);
    shp = cellfun(@(x) x.show, shp, 'UniformOutput', false);
    drawnow
    dt = toc; % get current timecost

    f.set('Name',['Time=', num2str(t,3), ' (s); ('...
      , 'FPS=',num2str(1/dt,3),')']); % 
end
