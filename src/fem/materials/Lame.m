function [Lambda,Mu] = Lame(E0,Nu0)
Lambda = (Nu0*E0)/((1+Nu0)*(1-2*Nu0));
Mu = E0/(2*(1+Nu0));
end