function y = polyq5(t, q0, qf, v0, vf, a0, af)
% CUBICPOLYTRAJECTORY Generates an output array of the coefficients ai
%  This functions solves for a quintic polynomial trajectory coefficents by taking in
%   start, t0, and end times, tf, and start, v0, end velocity ,vf, and
%   start q0 and end positions qf, start acceleration a0, end acceleration af

t0 = 0;
tf = 1;

A = [1 t0 (t0^2) (t0^3) (t0^4) (t0^5); 0 1 2*t0 3*(t0^2) 4*(t0^3) 5*(t0^4); ...
    0 0 2 6*t0 12*(t0^2) 20*(t0^3); 1 tf (tf^2) (tf^3) (tf^4) (tf^5); ...
    0 1 2*tf 3*(tf^2) 4*(tf^3) 5*(tf^4); 0 0 2 6*tf 12*(tf^2) 20*(tf^3)];

b = [q0; v0; a0; qf; vf; af];

x = A\b;

y = x(1) + x(2)*t + x(3)*t.^2 + x(4)*t.^3 + x(5)*t.^4 + x(6)*t.^5;

end
