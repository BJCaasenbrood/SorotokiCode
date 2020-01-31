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
        SDF;
        BdBox;
        
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
    
    obj.TextureStretch = 1;
    obj.Quality = 80;
    obj.FlipNormals = false;
    obj.AOBias = 0.01;
    obj.AOBit = 3;
    obj.AORadius = 0.3;
    obj.AOInvert = false;
            
    obj.AmbientOcclusion = false;
    obj.SubSurfaceScattering = false;
    obj.SSSPower = 2.2;
    obj.AOPower = 1;
    obj.Occlusion = [.2,.2,.2];
    
    for ii = 3:2:length(varargin)
        obj.(varargin{ii}) = varargin{ii+1};
    end

    
    obj = GenerateObject(obj,varargin);   
    
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
        'none','FaceVertexCData',Gmodel.TextureMap,'FaceColor','interp');
    
    set(gcf,'color',[255, 255, 255]/255); material dull;
    axis equal; axis(Gmodel.BdBox); axis off;
    daspect([1,1,1]);
    
    ax = gca; ax.Clipping = 'off';    
   
    Gmodel.FigHandle = hp;
    Gmodel.Figure = H;
    
    h.ActionPostCallback  = @myprecallback;
    Gmodel.Figure.UserData = Gmodel;
    
    function myprecallback(src,evnt)
        update(src.UserData);
    end

end

%--------------------------------------------------------------------- show
function update(Gmodel,varargin)
    
    Gmodel = updateNode(Gmodel);
    
    Gmodel = BakeCubemap(Gmodel,Gmodel.Texture);
    Gmodel.CameraPos = campos;
    
    Gmodel = updateTexture(Gmodel);
    
    set(Gmodel.FigHandle,'FaceVertexCData',Gmodel.TextureMap);
    drawnow;
end

%--------------------------------------------------------------------- show
function Gmodel = resetNode(Gmodel)
    Gmodel.Node = Gmodel.Node0;
end

%--------------------------------------------------------------------- show
function Gmodel = updateNode(Gmodel,varargin)
    
    if ~isempty(varargin), V = varargin{1};
    else, V = Gmodel.Node; end
    
    [Gmodel.VNormal,Gmodel.Normal] = TriangleNormal(V,...
        Gmodel.Element);

    set(Gmodel.FigHandle,'Vertices',Gmodel.Node);
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
    CC = ((1-Gmodel.AOTextureMap).^p);
    
    for ii = 1:length(Gmodel.Node)
        N(ii,:) = ColorMultiply(N(ii,:),repmat(CC(ii)+(1-CC(ii)),1,3));
        N(ii,:) = AffineColorMix(N(ii,:),E,...
            [0,0,0],[(CC(ii)),1-CC(ii),0]);
    end
    end
        
    if Gmodel.AmbientOcclusion
    E = Gmodel.Occlusion;
    p = Gmodel.AOPower;
    CC = ((Gmodel.AOTextureMap).^p);
    
    for ii = 1:length(Gmodel.Node)
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
    
    if Gmodel.SubSurfaceScattering, Gmodel.AOInvert = true; end
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
    elseif isa(msh{1},'function_handle')
       if length(msh{2}) ~= 6
           cout('err','* specify a correct bounding box of size 6 x 1');
       end
              
       Gmodel.SDF = msh{1};
       Gmodel.BdBox = msh{2};
       [X,Y,Z,P] = MeshGridding(Gmodel.BdBox,Gmodel.Quality);
       P = single(P);
       D = Gmodel.SDF(single(P));
       D = reshape(D(:,end),size(X));
       
       [f,v,~] = MarchingCubes(single(X),single(Y),...
           single(Z),single(D),-1e-6);
       
       %f = fliplr(f);
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
    Gmodel.BdBox = BoundingBox(Gmodel.Node); 
    Gmodel.TextureStretch = 1.0;
    
    fcell = num2cell(Gmodel.Element,2);
    [~,Gmodel.V2F,Gmodel.F2V] = ElementAdjecency(fcell);

end

%-------------------------------------------------- generate graphics model
function Gmodel = BakeCubemap(Gmodel,Cubemap)
      
    Ny = size(Cubemap,1); 
    Nx = size(Cubemap,2);
    
    if ishandle(101), Phi = view;
    else, Phi = eye(4); end
    
    Phi = transpose(Phi(1:3,1:3));
    
    if Gmodel.FlipNormals, Normals = -Gmodel.VNormal*Phi;
    else, Normals = -Gmodel.VNormal*Phi; end
 
    alpha = 1/Gmodel.TextureStretch;
    UV = SphereMapping(Normals,alpha);
       
    EnviromentReflect = zeros(length(Gmodel.Node),3);
    rgbImage = fliplr((Cubemap));
      
    R = rgbImage(:,:,1);
    G = rgbImage(:,:,2);
    B = rgbImage(:,:,3);

    y = (clamp(round(Nx*real(UV(:,1))),1,Nx));
    x = (clamp(round(Ny*real(UV(:,2))),1,Ny));
    
    ind = sub2ind([Nx Ny],y(:),x(:));
    EnviromentReflect(:,1) = double(R(ind))/255;
    EnviromentReflect(:,2) = double(G(ind))/255;
    EnviromentReflect(:,3) = double(B(ind))/255;
    
    RGB = zeros(length(Gmodel.Node),3);
%     
    for ii = 1:3
        RGB(:,ii) = (EnviromentReflect(:,ii));
    end
        
    Gmodel.TextureMap = RGB;    
end

%-------------------------------------------------- generate graphics model
function Gmodel = BakeAmbientOcclusion(Gmodel)

Bias = Gmodel.AOBias;
Res = round(2.^Gmodel.AOBit);
BB = BoundingBox(Gmodel.Node); 
R = Gmodel.AORadius*mean([abs(BB(2)-BB(1)),abs(BB(4)-BB(3)),...
    abs(BB(6)-BB(5))]);

Gmodel.BdBox = BB;

Normals = Gmodel.VNormal;
Centers = Gmodel.Node;
AO = zeros(length(Centers),1);
        
if Gmodel.AOInvert, dir = -1; else, dir = 1; end

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

f = fv.faces;
Gmodel.AOTextureMap = TextureSmoothing(f,AO,5);
end
    
end
end

