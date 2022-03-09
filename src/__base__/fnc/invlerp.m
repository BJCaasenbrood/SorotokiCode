function y = invlerp(a,b,v)
%LERP.m linear extrapolation between A and B: y = A + (B-A)x    x:=[0,1]
y = (v-a)/(b-a);
end
