function y = pwl(x,n,c)
% n = index
% c = number of cuts
if n > c || n < 1 
   error('Index number cannot be large that number of piece-wise segments') 
end
d = (n-1)/(c-1);    
s = (-abs((x-d)*(c-1))+1);

% s1 = (x >= (n-1)/c) & (x < (n)/c);
% s2 = (x == (n)/c);
% S = s1 | s2;

if c == 1
    y = double(x)*0 + 1;
else
    y = max(double(s),0);
end
end

