function E = Texture2Emission(C,factor)
if nargin < 2, factor = 1; end
Csum = cumsum(C,2);
[~,I] = max(Csum(:,end));

E = clamp(factor*C(I,:),0,1);
end

