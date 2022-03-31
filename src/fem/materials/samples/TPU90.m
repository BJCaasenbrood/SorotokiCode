function Material = TPU90(Nu)
if nargin< 1, Nu = 0.45; end
Material=  NeoHookeanMaterial('E',6.9,'Nu',Nu);

Material.Zeta = 0.1;
Material.Rho  = 1.07e-9;
end

