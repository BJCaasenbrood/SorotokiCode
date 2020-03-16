classdef Gmodel < handle

    properties (Access = public)
        Node;
        Element;
        NNode;
        NElem;
        Normal;
        VNormal;
        Texture;
        Emission;
        Occlusion;
    end
    
    properties (Access = private)
        Name;
        Node0;
        SDF;
        BdBox;
        
        Figure;
        FigHandle;
        TextureMap;
        TextureStretch;
        Shading;
        
        FlipNormals;
        AOPower;
        AOBias;
        AOBit;
        AORadius;
        AOInvert;
        AOTextureMap;
        
        SSSPower;
                
        AO;
        SSS;
        
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
    
    obj.Texture = grey;
    obj.TextureStretch = 1;
    obj.Quality = 80;
    obj.FlipNormals = false;
    obj.AOBias = 0.01;
    obj.AOBit = 3;
    obj.AORadius = 0.3;
    obj.AOInvert = false;
            
    obj.AO = false;
    obj.SSS = false;
    obj.SSSPower = 2.2;
    obj.AOPower = 1;
    obj.Occlusion = [.2,.2,.2];
    obj.Shading = 'Vertex';
    
    for ii = 3:2:length(varargin)
        obj.(varargin{ii}) = varargin{ii+1};
    end

    obj = GenerateObject(obj,varargin);   
    obj.NElem = length(obj.Element);
    obj.NNode = length(obj.Node);
    
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

%-------------------------------------------------------------- plot ground
function Gmodel = ground(Gmodel)
    Groundplane(Gmodel);
end

%--------------------------------------------------------------------- show
function Gmodel = render(Gmodel,varargin)
    
    H = figure(101);
    h = rotate3d;
    h.Enable = 'on';
    
    if strcmp(Gmodel.Shading,'Face'), shading = 'flat'; 
    else, shading = 'interp';
    end
    
    hp = patch('Vertices',Gmodel.Node,'Faces',Gmodel.Element,'linestyle',...
        'none','FaceVertexCData',Gmodel.TextureMap,'FaceColor',shading);
    
    set(gcf,'color',[255, 255, 255]/255); material dull;
    axis equal; axis(Gmodel.BdBox); axis off;
    daspect([1,1,1]);
    
    ax = gca; ax.Clipping = 'off';    
   
    Gmodel.FigHandle = hp;
    Gmodel.Figure = H;
    
    h.ActionPostCallback  = @myprecallback;
    Gmodel.Figure.UserData = Gmodel;
    
    %update(Gmodel);
    %drawnow;
    
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
      
    if strcmp(Gmodel.Shading,'Face'), shading = 'flat'; 
    else, shading = 'interp';
    end
    
    set(Gmodel.FigHandle,'FaceVertexCData',Gmodel.TextureMap,...
        'Facecolor',shading);
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
    
    if strcmp(Gmodel.Shading,'Face'), M = length(Gmodel.Element);
    else, M = length(Gmodel.Node);
    end
    N = Gmodel.TextureMap;
    
    if Gmodel.SSS
    if isempty(Gmodel.Emission)
        E = Texture2Emission(Gmodel.TextureMap,1.1);
        Gmodel.Emission = E;
    else, E = Gmodel.Emission;
    end
    
    p = Gmodel.SSSPower;
    CC = ((1-Gmodel.AOTextureMap).^p);
    
    for ii = 1:M
        N(ii,:) = ColorMultiply(N(ii,:),repmat(CC(ii)+(1-CC(ii)),1,3));
        N(ii,:) = AffineColorMix(N(ii,:),E,...
            [0,0,0],[(CC(ii)),1-CC(ii),0]);
    end
    end
        
    if Gmodel.AO
    E = Gmodel.Occlusion;
    p = Gmodel.AOPower;
    CC = ((Gmodel.AOTextureMap).^p);
    
    for ii = 1:M
        N(ii,:) = ColorMultiply(N(ii,:),repmat(CC(ii)+(1-CC(ii)),1,3));
        N(ii,:) = AffineColorMix(N(ii,:),E,...
            [0,0,0],[(CC(ii)),1-CC(ii),0]);
    end
    end
    
    if strcmp(Gmodel.Shading,'face')
       N =  TextureSmoothing(Gmodel.Element,N,5);
    end
    
    Gmodel.TextureMap = N;
    
end

%--------------------------------------------------------------------- show
function showAO(Gmodel,varargin)

    AOmap = TextureSmoothing(Gmodel.Element,Gmodel.AOTextureMap,10);
    
    set(Gmodel.FigHandle,'FaceVertexCData',AOmap,'facecolor','interp');
    colormap(gray);
    drawnow;
end

%--------------------------------------------------------------------- bake
function Gmodel = bake(Gmodel)
    
    Gmodel = BakeCubemap(Gmodel,Gmodel.Texture);
    
    if Gmodel.SSS, Gmodel.AOInvert = true; end
    if Gmodel.AO || Gmodel.SSS
    Gmodel = BakeAmbientOcclusion(Gmodel);
    end
end

%---------------------------------------------------------------- export
function export(Gmodel,type)
if nargin == 1, type = 'stl'; end

fv = struct;
f = Gmodel.Element;
v = Gmodel.Node;

[v, ~, indexn] =  unique(v, 'rows');
f = indexn(f);

if strcmp(type,'stl')
fv.vertices = v;
fv.faces = f;
stlwriter('test.stl',fv);
elseif strcmp(type,'obj')
fv.vertices = v;
fv.faces = f;
fv.objects(1).type='f';
fv.objects(1).data.vertices=fv.faces;
objwriter(fv,'test.obj');

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
       elseif strcmp(ext,'.obj'), [v,f] = objreader(msh{1});
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
          single(Z),single(D),1e-6);
%        [f,v] = isosurface(single(X),single(Y),...
%            single(Z),single(D),1e-6);
       
       %f = fliplr(f);    
    elseif length(msh) == 2
       v = msh{1};
       f = msh{2};        
    end
    
    [vn,fn] = TriangleNormal(v,f);
    
    Gmodel.Node = v;
    Gmodel.Node0 = v;
    Gmodel.Element = f;
    Gmodel.VNormal = -vn;
    Gmodel.Normal = fn;  
    Gmodel.RMatrix = eye(4);
    Gmodel.BdBox = BoundingBox(Gmodel.Node); 
    Gmodel.TextureStretch = 1.0;
    
    fcell = num2cell(Gmodel.Element,2);
    [~,Gmodel.V2F,Gmodel.F2V] = ElementAdjecency(fcell);
    
    Gmodel.TextureMap = Normal2RGB(Gmodel.Normal);

end

%-------------------------------------------------- generate graphics model
function Gmodel = BakeCubemap(Gmodel,Cubemap)
      
    Ny = size(Cubemap,1); 
    Nx = size(Cubemap,2);
    
    if ishandle(101), Phi = view;
    else, Phi = eye(4); end
    
    Phi = transpose(Phi(1:3,1:3));
    
    if strcmp(Gmodel.Shading,'Face')
        if Gmodel.FlipNormals, Normals = -Gmodel.Normal*Phi;
        else, Normals = -Gmodel.Normal*Phi; end
        N = length(Gmodel.Element);
    else
        if Gmodel.FlipNormals, Normals = -Gmodel.VNormal*Phi;
        else, Normals = -Gmodel.VNormal*Phi; end
        N = length(Gmodel.Node);
    end
 
    alpha = 1/Gmodel.TextureStretch;
    UV = SphereMapping(Normals,alpha);
       
    EnviromentReflect = zeros(N,3);
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
    
    RGB = zeros(N,3);
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

% if strcmp(Gmodel.Shading,'Face')
% Centers = zeros(length(Gmodel.Element),3);
% for el = 1:length(Gmodel.Element)
%     Centers(el,:) = mean(Gmodel.Node(Gmodel.Element(el,:),:),1);
% end
% Normals = Gmodel.Normal;
% else
Normals = Gmodel.VNormal;
Centers = Gmodel.Node;
% end

AOmap = zeros(length(Centers),1);
        
if Gmodel.AOInvert, dir = 1; else, dir = -1; end

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

[vnew, ~, indexn] =  unique(Gmodel.Node, 'rows');
fnew = indexn(Gmodel.Element);

fv = struct; 
fv.faces = fnew; 
fv.vertices =  vnew;

INC = inpolyhedron(fv,Cpts);
PolyFlag = INC(ic);
idx = 0;

for ii = 1:length(Centers)
    flag = PolyFlag(1+idx:Res+idx);
    if Gmodel.AOInvert, AOmap(ii) = 1 - sum(flag)/Res;
    else, AOmap(ii) = 1 - sum(flag)/Res;
    end
    
    idx = idx + Res;
end

f = fv.faces;
if ~strcmp(Gmodel.Shading,'Face')
Gmodel.AOTextureMap = TextureSmoothing(f,AOmap,5);
else
Gmodel.AOTextureMap = AOmap;

end
end

function Groundplane(Gmodel)
tmp = Gmodel.BdBox; a = 0.1;
tmp(1) = Gmodel.BdBox(1)-a*( Gmodel.BdBox(2) - Gmodel.BdBox(1)); 
tmp(2) = Gmodel.BdBox(2)+a*( Gmodel.BdBox(2) - Gmodel.BdBox(1)); 
tmp(3) = Gmodel.BdBox(3)-a*( Gmodel.BdBox(4) - Gmodel.BdBox(3)); 
tmp(4) = Gmodel.BdBox(4)+a*( Gmodel.BdBox(4) - Gmodel.BdBox(3)); 
tmp(5) = Gmodel.BdBox(5)-0*( Gmodel.BdBox(6) - Gmodel.BdBox(5)); 
tmp(6) = Gmodel.BdBox(6)+a*( Gmodel.BdBox(6) - Gmodel.BdBox(5)); 
B = tmp;

Nx = 4;
Ny = 4;

x = linspace(B(1),B(2),Nx+1);
y = linspace(B(3),B(4),Ny+1);

[X,Y] = meshgrid(x,y);

v = [tmp(1),tmp(3), tmp(5);  
     tmp(2),tmp(3), tmp(5);
     tmp(2),tmp(4), tmp(5);
     tmp(1),tmp(4), tmp(5)];
 
f = [1,2,3,4];

[I,~] = imread('checker.jpg');

dX = (tmp(2)-tmp(1));
dY = (tmp(4)-tmp(3));

if dY/dX >= 2, I = vertzcat(I,I);
elseif dX/dY >= 2, I = horzcat(I,I);
end

hold all
warpim(X,Y,X*0 + tmp(5),I);
hold off;

patch('Faces',f,'Vertices',v,...
    'Linewidth',1.5,'linestyle','-','FaceColor','none',...
    'EdgeColor',[1 1 1]*0.5);
end
    
end
end

