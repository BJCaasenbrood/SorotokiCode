clr; 
%% load hand preset
[shp, obj] = preset.shapes.softhand;
f = figure(101);

% dt = 1/30;
while true
    axis([-70,70,-70,70,-70,90]);

    t = shp{1}.solver.Time;
    for ii = 1:5, shp{ii}.system.fControl = [500; 0] * smoothstep(0.5 * t ); end

    tic;
    shp = cellfun(@(x) x.update(1/60), shp, 'UniformOutput', false);
    shp = cellfun(@showRenderShapes, shp, 'UniformOutput', false);
    dt = toc;

    disp(dt)
    drawnow limitrate;
    f.set('Name',['FPS=',num2str(1/dt,3)]);

    if norm(shp{1}.solver.sol.dx(:)) < 1e-6; break; end
end