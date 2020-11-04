function y = lerp(a,b,t)
%LERP.m linear extrapolation between A and B: y = A + (B-A)x    x:=[0,1]
if t > 1, t = 1;
elseif t < 0, t = 0;
end
y = a+t*(b-a);
end
