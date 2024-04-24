clr;
% generate fem model from diamond bot
msh = preset.mesh.diamond_bot;
fem = Fem(msh,'isNonlinear',true);
fem = fem.addSupport('bottom',[1,1]);

% assigning materials
mat = NeoHookean(2.5, 0.35);
fem = fem.addMaterial(mat);

fem.options.LineStyle = 'none';
fem = fem.eigen(); 
Ia = fem.system.Ia;

figure(101);
for ii = 1:6
    subplot(2,3,ii);

    % overwrite (exaggerated) solutions
    fem.solver.sol.x = fem.solver.sol.pod.V(:,ii) * (350 / sqrt(ii) );
    fem = fem.compute(); % compute internal stresses

    showVonMisesFem(fem);
    axis( [-20 185 0 120] );
    box on; axis on;
end
