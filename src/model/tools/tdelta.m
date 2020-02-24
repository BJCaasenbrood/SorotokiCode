function y = tdelta(x,h)

y = h^(-1)*max(1-abs(x/h),0);

end

