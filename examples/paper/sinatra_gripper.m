clr;
% get preset model
[shp, obj] = preset.shapes.sinatra_gripper;

% set axis for figure;
axis([-70,70,-70,70,-70,90]);

% simulation loop
while true

    % get time;
    t = shp{1}.solver.Time;

    % assgin acutation to each finger
    for ii = 1:6, shp{ii}.system.fControl = ...
        [0.5; 0; 0] * smoothstep(0.15 * t ); 
    end

    % update dynamics and render shapes
    shp = cellfun(@(x) x.update(), shp, 'UniformOutput', false);
    shp = cellfun(@showRenderShapes, shp, 'UniformOutput', false);

    % stop if zero-velocity is reached
    if norm(shp{1}.solver.sol.dx(:)) < 1e-6
        break; 
    end
end