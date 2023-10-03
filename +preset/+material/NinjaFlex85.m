function Material = NinjaFlex85(D)
%http://journals.ub.uni-magdeburg.de/ubjournals/index.php/techmech/article/view/564/540
%Material=  NeoHookeanMaterial('E',12,'Nu',.4);
%https://mds.marshall.edu/cgi/viewcontent.cgi?article=2251&context=etd
if nargin < 1, D = 10; end
Material = Mooney('C10',0.77,'C01',2.94);
end

