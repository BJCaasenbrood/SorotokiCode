function y = clamp(x,low,upp)
%CLAMP return bounded value between a lower-bound and upper-bound. The
%   function reads as y = clamp(x,low,upp) with low and upp the lower and
%   upper bound of x, respectively. Usage:
%
%   t = linspace(0,2*pi,500);
%   y = clamp(2*sin(t),-1,1);
%   plot(t,y);
%
%   Last edit: July 15, 2022.

y=min(max(x,low),upp);

end