function y = clamp(x,low,upp)
% function clamp(x,low,upp) - SOROTOKI
%-------------------------------------------------------------------------
% Description:
% Return bounded value between a lower-bound and upper-bound
%
% Usage:
%   t = linspace(0,2*pi,500);
%   y = clamp(2*sin(t),-1,1);
%   plot(t,y);
%--------------------------------------------------------------------------

y=min(max(x,low),upp);

end