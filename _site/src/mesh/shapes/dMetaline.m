function d = dMetaline(P,v1,v2,R,w)
d = zeros(length(P),1);

if length(R) == 1, R = d + R; end
if length(w) == 1, w = d + w; end

r = PointToLine(P, v1, v2);
d = SoftObject(r,w,R) + d; 

d = d+1e-3;
d = [d,d];
end

%-------------------------------------------------------------------------%
function d = RawGaussian(r,a,b)
 a = a(1); b = b(1);
 d = -b./r(:).^2;
 d = d - max(min(exp(1/a),5000),1)*mean(maxk(d,10));
end

%-------------------------------------------------------------------------%
function d = SoftObject(r,a,b)

 a = a(1); b = b(1);

 I1 = (r<=b);
 I2 = (r>b);
 
% d = zeros(length(r),1);
 %R = r(I1);
 R = r;
 d = -a*(1 - (4*R.^6)./(9*b^6)+(17*R.^4)./(9*b^4)-(22*R.^2)/(9*b^2));
 d(I2) = 1e-6*d(I2) + 1e-6;
 
% d = -1./r(:).^2
% ;
end


function d = PointToLine(pt, v1, v2)
% pt should be nx3
% v1 and v2 are vertices on the line (each 1x3)
% d is a nx1 vector with the orthogonal distances
V1 = repmat(v1,size(pt,1),1);
V2 = repmat(v2,size(pt,1),1);
a = V1 - V2;
b = pt - V2;
d = sqrt(sum(cross(a,b,2).^2,2))./sqrt(sum(a.^2,2));

for ii = 1:length(pt)
dse = sqrt(sum((v1 - v2).^2));
dsp = sqrt(sum((v1 - pt(ii,:)).^2));
dep = sqrt(sum((v2 - pt(ii,:)).^2));

if sqrt(dse^2 + dep^2) >= dsp, d(ii) = dep;
elseif sqrt(dse^2 + dsp^2) >= dep, d(ii) = dsp;
end
end

end