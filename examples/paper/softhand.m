clr; 
% load hand preset
[shp, obj] = preset.shapes.softhand;
f = figure(101);

dt = 1/40; 
t = 0;

% simulation loop
while t <= 11.0
    t = shp{1}.solver.Time;
    for ii = 1:5
      shp{ii}.system.fControl = [500; 0] * smoothstep(0.15 * t ); 
    end

    tic;
    shp = cellfun(@(x) x.update(dt), shp, 'UniformOutput', false);
    shp = cellfun(@(x) x.show, shp, 'UniformOutput', false);
    dt = toc;
end
