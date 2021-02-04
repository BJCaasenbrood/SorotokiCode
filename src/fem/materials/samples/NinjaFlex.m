function Material = NinjaFlex
%http://journals.ub.uni-magdeburg.de/ubjournals/index.php/techmech/article/view/564/540
%Material=  NeoHookeanMaterial('E',12,'Nu',.4);

Material = MooneyMaterial('C10',0.77,'C01',2.94);
end

