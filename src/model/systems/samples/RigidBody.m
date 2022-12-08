function sys = RigidBody(M,varargin)

q0 = [1,0,1,0];
s0 = [100,0,0];
v0 = [0,0,0];
w0 = [0,0,0];

sys = StateSpace({@(x,u,t) FLOW(x,u,t,M)},...
                 [q0,s0,w0,v0],'NInput',6);

end

function dx = FLOW(x,u,t,M)
q = x(1:4);
s = x(5:7);
w = x(8:10);
v = x(11:13);

I = M(4:6,4:6);
J = M(1:3,1:3);

Uw = u(1:3);
Uv = u(4:6);

Q = (q)/norm(q);
R = quat2rot(Q);

ds = v;                         % position update
dv = I\(R*Uv);                  % velocity update
dq = 0.5*OmegaOperator(w)*Q;    % quaternion update
dw = J\(Uw-skew(w)*(J*w));      % angular velocity update

dx = [dq;ds;dw;dv];
end
%--------------------------------------------------------------------------
function A = OmegaOperator(w)
A = zeros(4);
A(1,2:4)   = -w.';
A(2:4,1)   = +w;
A(2:4,2:4) = skew(w);
end
%--------------------------------------------------------------------------
function y = skew(v)
y = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0] ;
end