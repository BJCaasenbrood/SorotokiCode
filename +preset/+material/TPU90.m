function Material = TPU90(Nu)
if nargin< 1, Nu = 0.4; end
Material=  NeoHookeanMaterial('E',6.9,'Nu',Nu);

Material.params.Rho = 1070e-12;
Material.params.Zeta = 0.1;
end

