function y = pcc(x,n,c)
% n = index
% c = number of cuts
if n > c || n < 1 
   error('Index number cannot be large that number of piece-wise segments') 
end
    
s1 = (x >= (n-1)/c) & (x < (n)/c);
s2 = (x == (n)/c);
s = s1 | s2;

y = double(s);
end

