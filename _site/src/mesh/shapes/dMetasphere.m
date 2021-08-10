function d = dMetasphere(P,xc,yc,zc,R,w)

d = zeros(length(P),1);

% if length(R) == 1, R = d + R; end
% if length(w) == 1, w = d + w; end

for ii = 1:length(xc)
    r = sqrt((P(:,1)-xc(ii)).^2+(P(:,2)-yc(ii)).^2 + ...
        +(P(:,3)-zc(ii)).^2);
    d = RawGaussian(r,w(ii),R) + d; 
    %d = PiecewiseGaussian(r,w,R) + d; 
    %d = SoftObject(r,w,R) + d;
    %d = BlobbyMolecules(r,w,R) + d;
end

%d = d;
d = [d,d];

end

%-------------------------------------------------------------------------%
function d = RawGaussian(r,a,b)

 a = a(1); b = b(1);
 d = -b./r(:).^2;
 d = d - max(min(exp(1/a),5000),1)*mean(maxk(d,10));
end

%-------------------------------------------------------------------------%
function d = PiecewiseGaussian(r,a,b)

 a = a(1); b = b(1);

 I1 = (r<=b/3) & (r>=0);
 I2 = (r<=b) & (r>b/3);
 I3 = (r>b);
 
 d = zeros(length(r),1);
 d(I1) = -a*(1-3*(b.^-2)*r(I1).^2);
 d(I2) = -(3*a/2)*(1-r(I2)/b).^2;
 d(I3) = 1e-3;
 
 d = -1./r(:).^2;
 d = d - mean(maxk(d,10));
end

%-------------------------------------------------------------------------%
function d = BlobbyMolecules(r,a,b)

 a = a(1); 
 b = b(1);

 I1 = (r<=b);
 I2 = (r>b);
 
 d = zeros(length(r),1);
 expf = arrayfun(@(x) exp(-b*x^2),d);
 d = -a*expf;
% d = -1./r(:).^2
% ;
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