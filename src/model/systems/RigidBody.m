classdef RigidBody

    properties (Access = public)
        NDim;         % state dimensions
        NInput;       % input dimension
        Nonlinear;    % nonlinear true
        X0;           % initial conditions
       
        Gravity;      % gravity 
        
        Log;          % log file
        Gmodel;       % graphics model
        
    end
    
    properties (Access = private)
        Mtt;          % inertia tensor
    end
   
%--------------------------------------------------------------------------
methods  
%----------------------------------------------- MODAL SHAPE RECONSTRUCTION
function obj = RigidBody(varargin) 
    
    % construct linear ODE
    obj.Nonlinear = true;
    obj.NDim   = 13;
    obj.NInput = 6;
    obj.X0     = zeros(13,1);
    obj.X0(1)  = 1;
    obj.Gravity = zeros(3,1);
    
    for ii = 1:2:length(varargin)
        obj.(varargin{ii}) = varargin{ii+1};
    end
    
   % [~,obj] = flow(obj,obj.X0);
    
end
%---------------------------------------------------------------------- get     
function varargout = get(RigidBody,varargin)
    if nargin > 1
        varargout{nargin-1,1} = [];
        for ii = 1:length(varargin)
            varargout{ii,1} = Model.(varargin{ii});
        end
    else
        varargout = RigidBody.(varargin);
    end
end     
%--------------------------------------------------- get orientation matrix
function y = getMass(RigidBody)
    y = RigidBody.Mtt(4,4);
end
%--------------------------------------------------- get orientation matrix
function y = getOrientationMatrix(RigidBody)
    if isempty(RigidBody.Log)
        y = eye(3);
    else
        y = quat2rot(RigidBody.Log.X(end,1:4));
    end
end
%--------------------------------------------------- get orientation matrix
function y = getPosition(RigidBody)
    if isempty(RigidBody.Log)
        y = zeros(3,1);
    else
        y = RigidBody.Log.X(end,5:6).';
    end
end
%---------------------------------------------------------------------- set
function RigidBody = set(RigidBody,varargin)
    for ii = 1:2:length(varargin)
        RigidBody.(varargin{ii}) = varargin{ii+1};
    end
end
%---------------------------------------------------------------------- set
function RigidBody = setInitialSE3(RigidBody,x)
    RigidBody.X0(1:4) = rot2quat(x(1:3,1:3));
    RigidBody.X0(5:7) = x(1:3,4);
end
%---------------------------------------------------------------------- set
function RigidBody = setInitialVelocity(RigidBody,x)
    RigidBody.X0(11:13) = x;
end
%---------------------------------------------------------------------- set
function RigidBody = setInitialAngularVelocity(RigidBody,x)
    RigidBody.X0(8:10) = x;
end
%---------------------------------------------------------------------- set
function RigidBody = setGravity(RigidBody,varargin)
if isempty(varargin)
    varargin{1} = [0;0;-9810];
end

RigidBody.Gravity = varargin{1};
end
%------------------------------------------------------------------- render
function RigidBody = render(RigidBody,x)
    if isempty(RigidBody.Gmodel)
        RigidBody.Gmodel = Gmodel(sCube(1),'ShowProcess',false);
        RigidBody.Gmodel = Blender(RigidBody.Gmodel,'Translate',...
            [-0,-0,-1]);
        %RigidBody.Gmodel = Blender(RigidBody.Gmodel,'Scale',.5);
        RigidBody.Gmodel = Blender(RigidBody.Gmodel,'Fix');
        RigidBody.Gmodel.bake.render();
    end
    
    R = quat2rot(x(1:4));
    y = x(5:7);    
    RigidBody.Gmodel.reset();
    RigidBody.Gmodel = Blender(RigidBody.Gmodel,'SE3',SE3(R,y));
    RigidBody.Gmodel.update();
    
end
%-------------------------------------------------------------- return flow
function [dx, RigidBody] = flow(RigidBody,x,varargin)
    
    if ~isempty(varargin)
        u = varargin{1}(:);
        if numel(varargin)>1
            t = varargin{2};
        end
    else
        t = 0;
        u = 0;
    end
    
    [dx, RigidBody] = FLOW(RigidBody,x(:),u(:),t);
    
end
%----------------------------------------------------------- return hessian
function H = hessian(RigidBody,x,varargin)
    
   if ~isempty(varargin)
       u = varargin{1}(:);
       if numel(varargin)>1
           t = varargin{2};
       end
   else
       t = 0;
       u = 0;
   end 

   eps = 1e-6;
   H = zeros(RigidBody.NDim);
   f0 = flow(RigidBody,x,u,t);
   dx = eye(RigidBody.NDim)*eps;
   
   for ii = 1:RigidBody.NDim
       f1 = flow(RigidBody,x + dx(:,ii),u,t);
       H(:,ii) = (f1-f0)/eps;
   end

end

end
methods (Access = private)
%---------------------------------------------------------------------- set
function [dx, RigidBody] = FLOW(RigidBody,x,u,t)
q = x(1:4);
w = x(8:10);
v = x(11:13);

ag = RigidBody.Gravity(:);
Mt = RigidBody.Mtt;
I = Mt(4:6,4:6);
J = Mt(1:3,1:3);

Uw = u(1:3);
Uv = u(4:6);

Q = q;
Q = (q)/norm(q);                % ensure normalize
R = quat2rot(Q);

ds = v;                         % position update
dv = I\(R*Uv) + ag;             % velocity update
dq = 0.5*OmegaOperator(w)*Q;    % quaternion update
dw = J\(Uw-skew(w)*(J*w));      % angular velocity update

dx = [dq;ds;dw;dv];

% write log
RigidBody.Log.q   = Q;
RigidBody.Log.eta = [w(:);v(:)];
RigidBody.Log.g   = SE3(R,x(5:7));

    function y = skew(v)
        y = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0] ;
    end
    
    function A = OmegaOperator(w)
        A = zeros(4);
        A(1,2:4)   = -w.';
        A(2:4,1)   = +w;
        A(2:4,2:4) = skew(w);
    end
end

end
end

