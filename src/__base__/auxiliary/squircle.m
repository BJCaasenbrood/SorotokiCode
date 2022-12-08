function [x,y] = squircle(t,s,varargin)
%SQUIRCLE squircle(th,s,ra,rb)

% clamp s values
s = min(max(s,1e-6),1-1e-6);

if isempty(varargin)
    ra = 1; 
    rb = 1;
else
   if numel(varargin) == 1
       ra = varargin{1};
       rb = varargin{1};
   else
       ra = varargin{1};
       rb = varargin{2};
   end
end


RA = (ra)/(2*s);
RB = (rb)/(2*s);

Rx1 = RA*(2 + 2*s*sqrt(2) * cos(t) + s^2 * cos(2*t)).^0.5;
Ry1 = RB*(2 + 2*s*sqrt(2) * sin(t) - s^2 * cos(2*t)).^0.5;

Rx2 = -RA*(2 - 2*s*sqrt(2) * cos(t) + s^2 * cos(2*t)).^0.5;
Ry2 = -RB*(2 - 2*s*sqrt(2) * sin(t) - s^2 * cos(2*t)).^0.5;

x = Rx1 + Rx2;
y = Ry1 + Ry2;


end

