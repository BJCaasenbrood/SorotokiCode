function y = lerp(a,b,t)
%LERP.m linear extrapolation between A and B: y = A + (B-A)x    x:=[0,1]
y = a+t*(b-a);
end
