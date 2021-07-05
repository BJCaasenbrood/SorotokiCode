function Material = TPU90(Nu)
if nargin< 1, Nu = 0.33; end
Material=  NeoHookeanMaterial('E',69,'Nu',Nu);

Material.Zeta = 0.1;
Material.Rho  = 1.07e-9;
end

