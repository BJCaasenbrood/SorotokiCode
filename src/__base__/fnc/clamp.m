function y = clamp(x,low,upp)
% return bounded value clipped between bl and bu
y=min(max(x,low),upp);
end