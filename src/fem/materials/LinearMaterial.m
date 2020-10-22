classdef LinearMaterial

    properties (Access = public)
        Type = 'Linear';
        E = 1;
        Nu = 0.4;
        Mu;
        Lambda;
    end
   
%--------------------------------------------------------------------------
methods  
%--------------------------------------------------------------- Mesh Class
function obj = LinearMaterial(varargin) 
    
    for ii = 1:2:length(varargin)
        obj.(varargin{ii}) = varargin{ii+1};
    end
    
    obj = LameParameters(obj);
end

%---------------------------------------------------------------------- get     
function varargout = get(LinearMaterial,varargin)
    if nargin > 1
        varargout{nargin-1,1} = [];
        for ii = 1:length(varargin)
            varargout{ii,1} = LinearMaterial.(varargin{ii});
        end
    else
        varargout = LinearMaterial.(varargin);
    end
end
        
%---------------------------------------------------------------------- set
function LinearMaterial = set(LinearMaterial,varargin)
    for ii = 1:2:length(varargin)
        LinearMaterial.(varargin{ii}) = varargin{ii+1};
    end
end
    
%------------------------------ 2ND PIOLLA STRESSAND STIFFNESS FOR YEOH
function [S, D] = PiollaStress(LinearMaterial,F)
%Se = 2nd PK stress [S11, S22, S33, S12, S23, S13];
E0 = LinearMaterial.E;
Nu0 = LinearMaterial.Nu;
D=E0/(1-Nu0^2)*[1 Nu0 Nu0 0;Nu0 1 Nu0 0; Nu0 Nu0 1 0;0 0 0 (1-Nu0)/2];
D2=E0/(1-Nu0^2)*[1 Nu0 0; Nu0 1 0;0 0 (1-Nu0)/2];
C = F.'*F;

e = (1/2)*(C - eye(3));
S = D2*e;
end

%---------------------------------------------------------------------- set
function E = Emod(LinearMaterial)
   E = LinearMaterial.E;
end

end

methods (Access = private)
  
end
end

function LinearMaterial = LameParameters(LinearMaterial)
e = LinearMaterial.E;
v = LinearMaterial.Nu;

LinearMaterial.Lambda = (v*e)/((1+v)*(1-2*v));
LinearMaterial.Mu = e/(2*(1+v));
end


