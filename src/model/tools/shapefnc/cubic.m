function y = cubic(x,n,c)
% n = index
% c = number of cuts
if n > c || n < 1 
   error('Index number cannot be large that number of piece-wise segments') 
end
d = (n)/(c);    
s = (bstein((x-d)*(c-1),1,2));

% s1 = (x >= (n-1)/c) & (x < (n)/c);
% s2 = (x == (n)/c);
% S = s1 | s2;

if c == 1
    y = double(x)*0 + 1;
else
    y = max(double(s),0);
end
end


function y = bstein(x,v,n)
v = v-1;
v = n-1;
bnv = factorial(n)/(factorial(v)*factorial(n-v));
y = bnv*(x.^v).*(1-x).^(n-v);
end


