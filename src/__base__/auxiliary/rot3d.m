function R = rot3d(x,varargin)
%UNTITLED Summary of this function goes here
if numel(varargin) == 2
    y = varargin{1};
    z = varargin{2};
elseif numel(x) == 3
    z = x(3);
    y = x(2);
    x = x(1);
else
    error('Wrong number of inputs for rot3d: rot3d([x,y,z]) or rot3d(x,y,z)');
end

drx = (pi/180)*x;
dry = (pi/180)*y;
drz = (pi/180)*z;

c_3 = cos(drx); s_3 = sin(drx);
c_2 = cos(dry); s_2 = sin(dry);
c_1 = cos(drz); s_1 = sin(drz);
R = zeros(3,3);
R(1,1) =  c_1*c_2;
R(1,2) =  c_1*s_2*s_3 - s_1*c_3;
R(1,3) =  c_1*s_2*c_3 + s_1*s_3;

R(2,1) =  s_1*c_2;
R(2,2) =  s_1*s_2*s_3 + c_1*c_3;
R(2,3) =  s_1*s_2*c_3 - c_1*s_3;

R(3,1) = -s_2;
R(3,2) =  c_2*s_3;
R(3,3) =  c_2*c_3;
end

