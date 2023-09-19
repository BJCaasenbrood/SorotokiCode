clr; dt = 1/60;
%% load preset model
[shp, obj] = preset.shapes.suzumori_gripper();

% real-time simulatiom loop
while true

    % get time
    t = shp{1}.solver.Time;

    for ii = 1:4
        K = 4e-2 * shp{ii}.system.Stiffness;
        u = clamp(sign(sin(pi*t - (ii-1)*pi/8)), -Inf, 0);
        shp{ii}.system.fControl = K * [u; 0; 0; 0]; 
    end

    % update system and render
    tic;
    shp = cellfun(@(x) x.update(dt), shp, 'UniformOutput', false);
    shp = cellfun(@showRenderShapes, shp, 'UniformOutput', false);
    dt = toc;

    drawnow;
end