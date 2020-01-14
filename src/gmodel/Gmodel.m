classdef Gmodel < handle

    properties (Access = public)
        Name;
        Node;
        Node0;
        Element;
        Normal;
        VNormal;
        Texture;
        Emission;
        Occlusion;
    end
    
    properties (Access = private)
        Figure;
        FigHandle;
        TextureMap;
        TextureStretch;
        
        FlipNormals;
        AOPower;
        AOBias;
        AOBit;
        AORadius;
        AOInvert;
        AOTextureMap;
        
        SSSPower;
        
        AmbientOcclusion;
        SubSurfaceScattering;
        
        CameraPos;
        RMatrix;
        Quality;
        F2V;
        V2F;
    end
    
%--------------------------------------------------------------------------
methods  
%---------------------------------------------------------------- Fem Class
function obj = Gmodel(varargin) 
    
    obj = GenerateObject(obj,varargin);

    for ii = 3:2:length(varargin)
        obj.(varargin{ii}) = varargin{ii+1};
    end

    obj.TextureStretch = 1;
    obj.Quality = 20;
    obj.FlipNormals = false;
    obj.AOBias = 0.01;
    obj.AOBit = 3;
    obj.AORadius = 0.3;
    obj.AOInvert = true;
            
    obj.AmbientOcclusion = true;
    obj.SubSurfaceScattering = false;
    obj.SSSPower = 2.2;
    obj.AOPower = 1;
    obj.Occlusion = [.2,.2,.2];
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
function Gmodel = render(Gmodel,varargin)
    
    H = figure(101);
    h = rotate3d;
    h.Enable = 'on';
    hp = patch('Vertices',Gmodel.Node,'Faces',Gmodel.Element,'linestyle',...
        'none','FaceVertexCData',Gmodel.TextureMap,'FaceColor','flat');

    %hp.FaceColor = 'interp';
    set(gcf,'color',[255, 255, 255]/255);
    material dull;
    %camproj('perspective');
    axis equal;
    axis off;
    daspect([1,1,1]);
    
    ax = gca;               
    ax.Clipping = 'off';    
   
    Gmodel.FigHandle = hp;
    Gmodel.Figure = H;
    
    %set(h,'ButtonDownFilter',@myprecallback);
    h.ActionPostCallback  = @myprecallback;
    %h.ButtonDownFilter = @mypostcallback;
    Gmodel.Figure.UserData = Gmodel;
    
    function myprecallback(src,evnt)
        update(src.UserData);
    end

end

%--------------------------------------------------------------------- show
function update(Gmodel,varargin)
    
    %if sum(Gmodel.CameraPos) ~= sum(campos)
    Gmodel = BakeCubemap(Gmodel,Gmodel.Texture);
    Gmodel.CameraPos = campos;
    %end
    Gmodel = updateTexture(Gmodel);
    
    set(Gmodel.FigHandle,'FaceVertexCData',Gmodel.TextureMap);
    drawnow;
end

%--------------------------------------------------------------------- show
function Gmodel = resetNode(Gmodel)
    
    %if sum(Gmodel.CameraPos) ~= sum(campos)
    Gmodel.Node = Gmodel.Node0;
 
end

%--------------------------------------------------------------------- show
function updateNode(Gmodel,varargin)
    
    [vn,fn] = TriangleNormal(varargin{1},Gmodel.Element);

    Gmodel.VNormal = vn;
    Gmodel.Normal = fn;  
        
    set(Gmodel.FigHandle,'Vertices',varargin{1});
    %update(Gmodel);
end

%--------------------------------------------------------------------- show
function Gmodel = updateTexture(Gmodel)
    N = Gmodel.TextureMap;
    
    if Gmodel.SubSurfaceScattering
    if isempty(Gmodel.Emission)
        E = Texture2Emission(Gmodel.TextureMap,1.1);
        Gmodel.Emission = E;
    else, E = Gmodel.Emission;
    end
    
    p = Gmodel.SSSPower;
    CC = ((Gmodel.AOTextureMap).^p);
    
    for ii = 1:length(Gmodel.Element)
        N(ii,:) = ColorMultiply(N(ii,:),repmat(CC(ii)+(1-CC(ii)),1,3));
        N(ii,:) = AffineColorMix(N(ii,:),E,...
            [0,0,0],[(CC(ii)),1-CC(ii),0]);
    end
    end
        
    if Gmodel.AmbientOcclusion
    E = Gmodel.Occlusion;
    p = Gmodel.AOPower;
    CC = ((Gmodel.AOTextureMap).^p);
    
    for ii = 1:length(Gmodel.Element)
        N(ii,:) = ColorMultiply(N(ii,:),repmat(CC(ii)+(1-CC(ii)),1,3));
        N(ii,:) = AffineColorMix(N(ii,:),E,...
            [0,0,0],[(CC(ii)),1-CC(ii),0]);
    end
    end
    
    Gmodel.TextureMap = N;
end

%--------------------------------------------------------------------- show
function showAO(Gmodel,varargin)

    set(Gmodel.FigHandle,'FaceVertexCData',Gmodel.AOTextureMap);
    colormap(gray);
    drawnow;
end

%--------------------------------------------------------------------- bake
function Gmodel = bake(Gmodel)
    
    Gmodel = BakeCubemap(Gmodel,Gmodel.Texture);
    
    if Gmodel.AmbientOcclusion || Gmodel.SubSurfaceScattering
    Gmodel = BakeAmbientOcclusion(Gmodel);
    end
end
end


methods (Access = private)
    
%-------------------------------------------------- generate graphics model
function Gmodel = GenerateObject(Gmodel,varargin)
    
    msh = varargin{1};

    if isa(msh{1},'char')
       [~,Gmodel.Name,ext] = fileparts(msh{1});
       if strcmp(ext,'.stl'), [f,v] = stlreader(msh{1});
       elseif strcmp(ext,'.obj'), [f,v] = objreader(msh{1});
       else, cout('err','* extension not recognized');
       end
    elseif length(msh) == 2
       v = msh{1};
       f = msh{2};
    end
    
    [vn,fn] = TriangleNormal(v,f);
    
    Gmodel.Node = v;
    Gmodel.Node0 = v;
    Gmodel.Element = f;
    Gmodel.VNormal = vn;
    Gmodel.Normal = fn;  
    Gmodel.TextureMap = Normal2RGB(Gmodel.Normal);
    Gmodel.RMatrix = eye(4);
    
    fcell = num2cell(Gmodel.Element,2);
    [~,Gmodel.V2F,Gmodel.F2V] = ElementAdjecency(fcell);

end

%-------------------------------------------------- generate graphics model
function Gmodel = BakeCubemap(Gmodel,Cubemap)
      
    Ny = size(Cubemap,1); 
    Nx = size(Cubemap,2);
    
    Phi = view;
    Phi = transpose(Phi(1:3,1:3));
    
    if Gmodel.FlipNormals, Normals = -Gmodel.VNormal*Phi;
    else, Normals = -Gmodel.VNormal*Phi; end
 
    alpha = 1/Gmodel.TextureStretch;
    UV = SphereMapping(Normals,alpha);
       
    EnviromentReflect = zeros(length(Gmodel.Node),3);
    rgbImage = ((Cubemap));
      
    R = rgbImage(:,:,1);
    G = rgbImage(:,:,2);
    B = rgbImage(:,:,3);

    y = (clamp(round(Nx*real(UV(:,1))),1,Nx));
    x = (clamp(round(Ny*real(UV(:,2))),1,Ny));
    
%     for ii = 1:length(Gmodel.Node)
%         EnviromentReflect(ii,1) = double(rgbImage(y(ii), x(ii), 1))/255;
%         EnviromentReflect(ii,2) = double(rgbImage(y(ii), x(ii), 2))/255;
%         EnviromentReflect(ii,3) = double(rgbImage(y(ii), x(ii), 3))/255;
%     end
    ind = sub2ind([Nx Ny],y(:),x(:));
    EnviromentReflect(:,1) = double(R(ind))/255;
    EnviromentReflect(:,2) = double(G(ind))/255;
    EnviromentReflect(:,3) = double(B(ind))/255;
    %EnviromentReflect(I,2) = double(rgbImage(y(I), x(I), 2))/255;
    %EnviromentReflect(I,3) = double(rgbImage(y(I), x(I), 3))/255;
    
    RGB = zeros(length(Gmodel.Element),3);
%     
    for ii = 1:3
        RGB(:,ii) = Gmodel.V2F*(EnviromentReflect(:,ii));
    end
        
    Gmodel.TextureMap = RGB;    
%    Gmodel.TextureMap = EnviromentReflect;    
end

%-------------------------------------------------- generate graphics model
function Gmodel = BakeAmbientOcclusion(Gmodel)

Bias = Gmodel.AOBias;
Res = round(2.^Gmodel.AOBit);
BdBox = BoundingBox(Gmodel.Node); 
R = Gmodel.AORadius*mean([abs(BdBox(2)-BdBox(1)),abs(BdBox(4)-BdBox(3)),...
    abs(BdBox(6)-BdBox(5))]);

Normals = Gmodel.VNormal;
Centers = Gmodel.Node;
AO = zeros(length(Centers),1);
        
if ~Gmodel.AOInvert, dir = -1; else, dir = 1; end

npts = Res*length(Centers);
pts = zeros(npts,3);
idx = 0;

for ii = 1:length(Centers)
    T = GenerateHemiSphereBase(Res,Bias); 
    Rot = PlanarProjection(dir*Normals(ii,:));
    pts(1+idx:Res+idx,:) = R*((Rot*T.').') + Centers(ii,:);
    idx = idx + Res;
end

[Cpts,~,ic] = uniquetol(pts,R*0.5e-3,'ByRows',true);

%Cpts = pts;

fv = struct; 
fv.faces =  Gmodel.Element; 
fv.vertices =  Gmodel.Node;

INC = inpolyhedron(fv,Cpts);
PolyFlag = INC(ic);
idx = 0;

for ii = 1:length(Centers)
    flag = PolyFlag(1+idx:Res+idx);
    if Gmodel.AOInvert, AO(ii) = 1 - sum(flag)/Res;
    else, AO(ii) = sum(flag)/Res;
    end
    
    idx = idx + Res;
end

%AO = MeshSmoothing(Mesh,AO,5);
f = fv.faces;
Gmodel.AOTextureMap = Gmodel.V2F*TextureSmoothing(f,AO,5);
end
    
end
end

