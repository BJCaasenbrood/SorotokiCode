classdef Gmodel < handle

    properties (Access = public)
        Name;
        Node;
        Element;
        Normal;
        Texture;
    end
    
    properties (Access = private)
        TextureMap;
        CameraPos;
        Quality = 20;
        F2V;
        V2F;
    end
    
%--------------------------------------------------------------------------
methods  
%---------------------------------------------------------------- Fem Class
function obj = Gmodel(Input,varargin) 
    obj = GenerateObject(obj,Input);
    for ii = 1:2:length(varargin)
        obj.(varargin{ii}) = varargin{ii+1};
    end
    [~,obj.V2F,obj.F2V] = ElementAdjecency(num2cell(obj.Element,2));
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
    
    if sum(Gmodel.CameraPos) ~= sum(campos), Gmodel = bake(Gmodel); end
    
    cla;
    
    h = patch('Vertices',Gmodel.Node,'Faces',Gmodel.Element,'linestyle',...
        'none','FaceVertexCData',Gmodel.TextureMap,'FaceColor','flat');
    axis equal;
    
    %set(gcf,'color',[255, 255, 255]/255);
    material dull;
    %camproj('persp');
    axis off;
    drawnow;
end

%--------------------------------------------------------------------- bake
function Gmodel = bake(Gmodel)
    Gmodel = BakeCubemap(Gmodel,Gmodel.Texture);
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
    Gmodel.TextureMap = Normal2RGB(Gmodel.Normal);
end

%-------------------------------------------------- generate graphics model
function Gmodel = BakeCubemap(Gmodel,Cubemap)
      
    Ny = size(Cubemap,1); 
    Nx = size(Cubemap,2);
    
    Normals = VertexNormal(Gmodel.Node,Gmodel.Element);
    Normals = num2cell(Normals,2);
    Phi = view;
    Phi = Phi(1:3,1:3);
    
    UVMap = cellfun(@(x) SphereMapping(x,Phi),Normals,...
        'UniformOutput',false);
    
    EnviromentReflect = zeros(length(Gmodel.Node),3);
    rgbImage = (fliplr(Cubemap));
    
    for ii = 1:length(Gmodel.Node)
        
        Y = real(UVMap{ii}(1));
        X = real(UVMap{ii}(2));
        
        y = real(round(Nx*Y));
        x = real(round(Ny*X));
        [x,y] = CheckRGBRadial(x,y,Nx,Ny);
        
        EnviromentReflect(ii,1) = double(rgbImage(y, x, 1))/255;
        EnviromentReflect(ii,2) = double(rgbImage(y, x, 2))/255;
        EnviromentReflect(ii,3) = double(rgbImage(y, x, 3))/255;
    end
    
    RGB = zeros(length(Gmodel.Element),3);
%     
    for ii = 1:3
        RGB(:,ii) = Gmodel.V2F*(EnviromentReflect(:,ii));
    end
    
    Gmodel.TextureMap = RGB;
end
    
end
end

%------------------------------------------------------ SPHERICAL REMAPPING
function UV = SphereMapping(r,Rot)
r = r*Rot';
rx = r(1); ry = r(2); rz = r(3);%
m = 2*sqrt((rx)^2 + ry^2 + (abs(rz) + 1)^2);
v = 1.5*[(rx)/m, (ry)/m]; 

x = v(1); y = v(2);
UV = [y+0.5,x+0.5];
end

%------------------------------------------------- CHECK UV-TEXTURE MAPPING
function [x,y] = CheckRGBRadial(x,y,Nx,Ny)
R = floor(mean([Nx,Ny])/2);
Xm = floor(Nx/2); Ym = floor(Ny/2);
Delta = [x - Xm, y - Ym]; 
Theta = atan2(Delta(1),Delta(2));
dR = round(sqrt(Delta(1)^2 +  Delta(2)^2));

if dR >= R
   x = floor(Xm + sin(Theta) *0.99*R);
   y = floor(Ym + cos(Theta) *0.99*R);
end

if real(x) >= Nx, x = Nx-2; end
if real(y) >= Ny, y = Ny-2; end
if real(x) <= 1, x = 2; end
if real(y) <= 1, y = 2; end

if isnan(x), x = round(0.5*Nx)-5; end
if isnan(y), y = round(0.5*Ny)-5; end

x = real(x); y = real(y);
end

