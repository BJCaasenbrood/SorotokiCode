classdef Gmodel < handle

    properties (Access = public)
        Name;
        Node;
        Element;
        Normal;
        Texture;
    end
    
    properties (Access = private)
        Quality = 20;
    end
    
%--------------------------------------------------------------------------
methods  
%---------------------------------------------------------------- Fem Class
function obj = Gmodel(Input,varargin) 
    obj = GenerateObject(obj,Input);
    for ii = 1:2:length(varargin)
        obj.(varargin{ii}) = varargin{ii+1};
    end
end

%---------------------------------------------------------------------- get     
function varargout = get(Gmodel,varargin)
    if nargin > 1
        varargout{nargin-1,1} = [];
        for ii = 1:length(varargin)
            varargout{ii,1} = Gmodel.(varargin{ii});
        end
    else
        varargout = Gmodel.(varargin);
    end
end
        
%---------------------------------------------------------------------- set
function Gmodel = set(Gmodel,varargin)
    for ii = 1:2:length(varargin)
        Gmodel.(varargin{ii}) = varargin{ii+1};
    end
end

%--------------------------------------------------------------------- show
function h = show(Gmodel,varargin)
    h = patch('Vertices',Gmodel.Node,'Faces',Gmodel.Element,'linestyle',...
        'none','FaceVertexCData',Gmodel.Texture,'FaceColor','flat');
    axis equal;
    set(gcf,'color',[255, 255, 255]/255);
    material dull;
    camproj('persp');
    axis off;
    view(30,30);
    drawnow;
end

end
methods (Access = private)
    
%-------------------------------------------------- generate graphics model
function Gmodel = GenerateObject(Gmodel,Input)
    if isa(Input,'char')
       [~,Gmodel.Name,ext] = fileparts(Input);
       if strcmp(ext,'.stl'), [f,v] = stlreader(Input);
       elseif strcmp(ext,'.obj'), [f,v] = objreader(Input);
       else, cout('err','* extension not recognized');
       end
       
    end
    
    Gmodel.Node = v;
    Gmodel.Element = f;
    Gmodel.Normal = TriangleNormal(v,f);   
    Gmodel.Texture = Normal2RGB(Gmodel.Normal);
end
    
end
end


