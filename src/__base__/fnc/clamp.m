function y = clamp(x,low,upp)
% function clamp(x,low,upp)
%-------------------------------------------------------------------------
% Return bounded value clipped between lower-bound and upper-bound
%
% Usage:
%   t = linspace(0,2*pi,100);
%   y = clamp(2*sin(t),-2,3);
%--------------------------------------------------------------------------

y=min(max(x,low),upp);

end