classdef MDE
%Matrix Differential Equation Class for Cosserat Beams
    
    properties
        p;
        Phi;
        J;
    end
%--------------------------------------------------------------------------    
methods
%---------------------------------------------------------------------- set
function obj = MDE(n)
    obj.p   = zeros(3,1);
    obj.Phi = eye(3);
    obj.J   = zeros(6,n);
end
%---------------------------------------------------------------------- set
function Shapes = set(Shapes,varargin)
    for ii = 1:2:length(varargin)
        Shapes.(varargin{ii}) = varargin{ii+1};
    end
end
%------------------------------------------------------------ plus operator
function MDE = plus(MDE,mde)
    MDE.p   = MDE.p + mde.p;      % add the position vectors
    MDE.Phi = MDE.Phi + mde.Phi;  % add the orientation matrix
    MDE.J   = MDE.J + mde.J;      % add the Jacobian matrix
end
%----------------------------------------------------------- times operator        
function y = mtimes(a,b)
    if isa(b,'MDE')
        y     = b;
        y.p   = a*y.p;            % position vector times scalar
        y.Phi = a*y.Phi;          % orientation matrix times scalar
        y.J   = a*y.J ;           % Jacobian matrix times scalar
    else
        y     = a;
        y.p   = b*y.p;            % position vector times scalar
        y.Phi = b*y.Phi;          % orientation matrix times scalar
        y.J   = b*y.J ;           % Jacobian matrix times scalar
    end
end
%----------------------------------------------------------- times operator        
        
end
end

