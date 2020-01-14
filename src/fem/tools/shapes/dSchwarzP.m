function d = dSchwarzP(P,t,f)
if nargin < 3, f = 1/(2*pi); end
x = (f/2*pi + 1e-1)*P(:,1);
y = (f/2*pi + 1e-1)*P(:,2);
z = (f/2*pi + 1e-1)*P(:,3);
d = cos(x) + cos(y) + cos(z);
d = [d,abs(d)-t];
