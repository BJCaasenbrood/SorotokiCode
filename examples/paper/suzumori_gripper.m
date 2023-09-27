clr; dt = 1/220;
% load preset model

sdf = sCylinder(20,[-10,10],[0,0,-60]);
[shp, obj] = preset.shapes.suzumori_gripper('sdf',sdf);

% real-time simulatiom loop
while true

    % get time
    t = shp{1}.solver.Time;

    for ii = 1:4
        u = 2e-2 * smoothstep(t);
        K = shp{ii}.system.Stiffness;
        shp{ii}.system.fControl = K * [u; 0; 0; 0]; 
    end

    % update system and render
    tic;
    shp = cellfun(@(x) x.update(dt/2), shp, 'UniformOutput', false);
    shp = cellfun(@showRenderShapes, shp, 'UniformOutput', false);
    dt = toc;
    drawnow;
end