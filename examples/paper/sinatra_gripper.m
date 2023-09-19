
clr;
[shp, obj] = preset.shapes.sinatra_gripper;

while true
    axis([-70,70,-70,70,-70,90]);

    t = shp{1}.solver.Time;
    for ii = 1:6, shp{ii}.system.fControl = [0.5; 0; 0] * smoothstep(0.15 * t ); end

    shp = cellfun(@(x) x.update(), shp, 'UniformOutput', false);
    shp = cellfun(@showRenderShapes, shp, 'UniformOutput', false);

    if norm(shp{1}.solver.sol.dx(:)) < 1e-6; break; end
end