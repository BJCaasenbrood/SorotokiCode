clr;
% generate fem model from diamond bot
msh = preset.mesh.diamond_bot;
fem = Fem(msh,'isNonlinear',false);
fem = fem.addSupport('bottom',[1,1]);

% assigning materials
mat = NeoHookean(2.5, 0.35);
fem = fem.addMaterial(mat);

fem.options.LineStyle = 'none';
fem = fem.eigen(); 
Ia = fem.system.Ia;

figure(101);
for ii = 1:3
    subplot(1,3,ii);
    fem.solver.sol.x = fem.solver.pod.V(:,ii) * 250;
    fem.show('Linestyle','none');
    axis( [-20 185 0 120] );
end
