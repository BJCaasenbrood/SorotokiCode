function y = zclamp(x,low,upp)
% function zclamp(x,low,upp)
%-------------------------------------------------------------------------
% Return zero value clipped between lower-bound and upper-bound
%
% Usage:
%   t = linspace(0,2*pi,500);
%   y = zclamp(2*sin(t),-1,1);
%   plot(t,y);
%--------------------------------------------------------------------------
y = x;
for ii = 1: numel(x)
    if x(ii) < low
        y(ii) = 0;
    elseif x(ii) > upp
        y(ii) = 0;
    else
        y(ii) = x(ii);
    end
end

end