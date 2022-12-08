function y = bstein(x,v,n)
n = n-1;
bnv = factorial(n)/(factorial(v)*factorial(n-v));
y = bnv*(x.^v).*(1-x).^(n-v);
end

