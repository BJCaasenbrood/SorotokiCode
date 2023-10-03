% Credits to Matheus S.Xavier, Charbel D.Tawk, Yuen K.Yonga, Andrew J.Fleminga
% https://doi.org/10.1016/j.sna.2021.113199
function Material = Elastic50A(D)
if nargin < 1, D = .2; end
%  Material = YeohMaterial('C1',0.5150,'C2',0.0,'C3',0.0025,...
%      'D1',D,'D2',D,'D3',D);
Material = Yeoh('C1',0.5056,'C2',0.00,'C3',0.002562,...
     'D1',D,'D2',D,'D3',D);

Material.params.Rho  = 1200e-12;
Material.params.Zeta = 0.03;
Material.contact.NormalReaction = 0.3;
end

%C1 = 0.5243, C2 = 0, C3 = 0.004186
%C1 = 0.5056, C2 = 0, C3 = 0.002562