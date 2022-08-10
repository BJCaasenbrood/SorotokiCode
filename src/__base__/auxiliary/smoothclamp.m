function y = smoothclamp(x,low,upp,varargin)
%CLAMP return bounded value between a lower-bound and upper-bound. The
%   function reads as y = clamp(x,low,upp) with low and upp the lower and
%   upper bound of x, respectively. Usage:
%
%   t = linspace(0,2*pi,500);
%   y = clamp(2*sin(t),-1,1);
%   plot(t,y);
%
%   Last edit: July 15, 2022.

%y=min(max(x,low),upp);
if ~isempty(varargin)   
    lambda = varargin{1};
else
    lambda = 1;
end

t = lambda*(x-low)/(upp-low);
L1 = t < 0;
L2 = (t <= 1) & (~L1);
L3 = (~L1&~L2);

Y = (0)*L1 + (3*t.^2-2*t.^3).*L2 + 1*L3;

y = low + (upp-low)*Y;

% def smoothclamp(x, mi, mx): 
% return mi + (mx-mi)*(lambda t: 
% np.where(t < 0 , 0, np.where( t <= 1 , 3*t**2-2*t**3, 1 ) ) )(  )



end