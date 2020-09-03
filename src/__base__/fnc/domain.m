function D = domain(B0,B1,varargin)
%DOMAIN.m returns a bounded 2D/3D domain given input bounds

D = [B0,B1,B0,B1];
if nargin > 2
if varargin{1} == 3, D = [B0,B1,B0,B1,B0,B1];
end
end
end

