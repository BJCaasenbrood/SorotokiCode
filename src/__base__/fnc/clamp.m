function y = clamp(x,low,upp)
%CLAMP.m return bounded value clipped between lower bound and upper bound
y=min(max(x,low),upp);
end