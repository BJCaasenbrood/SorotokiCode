function y = tdelta(x,h)
y = h*max(1-abs(x*h),0);
end

