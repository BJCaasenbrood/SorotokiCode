classdef Gmodel < handle

    properties (Access = public)
        Node;
        Element;
        NNode;
        NElem;
        Texture;
        Emission;
        Occlusion;
        Alpha;
        View;
    end
    
    properties (Access = private)
        Name;
        Node0;
        Element0;
        SDF;
        BdBox;
        Center;
        
        Figure;
        FigHandle;
        FigAxis;
        TextureMap;
        TextureStretch;
        Shading;
        Colormap;
        LineStyle;
        LineColor;
        
        Normal;
        VNormal;
        FlipNormals;
        ShowAOActive;
        AOPower;
        AOBias;
        AOBit;
        AORadius;
        Radius;
        AOInvert;
        AOTextureMap;
        SSSTextureMap;
        AORenderComplete;
        SSSRenderComplete;
        
        SSSPower;
        SSSRadius;
                
        AO;
        SSS;
        
        Slice;
        
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
    
    obj.Texture = base;
    obj.TextureStretch = .85;
    obj.Quality = 80;
    obj.FlipNormals = false;
    obj.AOBias = 0.01;
    obj.AOBit = 3;
    obj.AORadius = 0.3;
    obj.SSSRadius = 0.3;
    obj.AOInvert = false;
    obj.Colormap = turbo;
    obj.LineStyle = 'none';
    obj.LineColor = col(1);
    obj.Alpha = 1.0;
            
    obj.AO = false;
    obj.SSS = false;
    obj.SSSPower = 2.2;
    obj.AOPower = 1;
    obj.Occlusion = [.2,.2,.2];
    obj.SSSRenderComplete = true;
    obj.AORenderComplete = true;
    obj.AOTextureMap = 1;
    obj.SSSTextureMap = 1;
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

%--------------------------------------------------------------------- copy
function obj = copy(Gmodel,varargin)
    if nargin > 2
       obj = Blender(Gmodel,varargin{:});
    else
       obj = Gmodel;
    end
    
end

%--------------------------------------------------------------------- copy
function Gmodel = fix(Gmodel)
        Gmodel.Node0 = Gmodel.Node;
end

%-------------------------------------------------------------- plot ground
function Gmodel = ground(Gmodel,gnd)
    if nargin < 2, Groundplane(Gmodel);
    else, Groundplane(Gmodel,gnd); end
end

%-------------------------------------------------------------- plot ground
function Gmodel = box(Gmodel)
    BoundingBox(Gmodel);
end

%--------------------------------------------------------------------- show
function Gmodel = render(Gmodel,varargin)
    
    if nargin < 2, H = figure(101);
    else, H = figure(varargin{1});
    end
    h = rotate3d;
    h.Enable = 'on';
    
%     if strcmp(Gmodel.Shading,'Face'), shd = 'flat'; 
%     else, shd = 'interp';
%     end
    
    hp = patch('Vertices',Gmodel.Node,'Faces',Gmodel.Element,'linestyle',...
        Gmodel.LineStyle,'edgecolor',Gmodel.LineColor,'FaceVertexCData',...
        Gmodel.TextureMap,'FaceColor','flat','FaceAlpha',Gmodel.Alpha);
    
    set(gcf,'color',gitpage); material dull;
    axis equal; axis(Gmodel.BdBox); axis off; %view(30,15);
    daspect([1,1,1]);
    
    ax = gca; 
    ax.Clipping = 'off';    
   
    Gmodel.FigHandle = hp;
    Gmodel.FigAxis = ax;
    Gmodel.Figure = H;
    
    h.ActionPostCallback  = @myprecallback;
    Gmodel.Figure.UserData = Gmodel;
    
    update(Gmodel);
    %drawnow;
    
    function myprecallback(src,evnt)
        class = whoClasses('Gmodel');
        for i = 1:length(class)
        	r = update(class{i},'tex');
        end
          class = whoClasses('Rig');
        for i = 1:length(class)
        	r = class{i}.update();
        end
    end

end

%--------------------------------------------------------------------- show
function vargout = update(Gmodel,varargin)
    
    if nargout < 1 && isempty(Gmodel.Slice)
        Gmodel = updateNode(Gmodel);
    end
    
    Gmodel = BakeCubemap(Gmodel,Gmodel.Texture);
    Gmodel.CameraPos = campos;
    
    Gmodel = updateTexture(Gmodel);
      
    if strcmp(Gmodel.Shading,'Face'), shading = 'flat'; 
    else, shading = 'interp';
    end
    
    if nargout < 1 || strcmp(varargin ,'tex')
    set(Gmodel.FigHandle,'FaceVertexCData',Gmodel.TextureMap,...
        'Facecolor',shading);
    
    vargout = [];
    else
       vargout{1} = Gmodel.FigHandle;
       vargout{2} = Gmodel.TextureMap;
       vargout{3} = shading;
       vargout{4} = Gmodel.Node;
    end
    
end

%--------------------------------------------------------------------- show
function Gmodel = reset(Gmodel)
    Gmodel.Node = Gmodel.Node0;
    Gmodel.Element = Gmodel.Element0;
end

%--------------------------------------------------------------------- show
function Gmodel = center(Gmodel)
    Gmodel.BdBox = boxhull(Gmodel.Node); 
end

%--------------------------------------------------------------------- show
function Gmodel = updateNode(Gmodel,varargin)
    
    if ~isempty(varargin), V = varargin{1};
    else, V = Gmodel.Node; end
    
    [Gmodel.VNormal,Gmodel.Normal] = TriangleNormal(V,...
        Gmodel.Element);

    set(Gmodel.FigHandle,'Vertices',Gmodel.Node);
    
    %drawnow limitrate;
end

%--------------------------------------------------------------------- show
function Gmodel = updateElements(Gmodel,varargin)
    set(Gmodel.FigHandle,'Faces',Gmodel.Element);
end

%--------------------------------------------------------------------- show
function [Gmodel,map] = updateTexture(Gmodel)
    
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
        CC = ((1-Gmodel.SSSTextureMap).^p);
        
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
    map = N;
    
end

%--------------------------------------------------------------------- show
function showMap(Gmodel,Request)
    
    switch(Request)
    case('AO');  P = TextureSmoothing(Gmodel.Element,1./Gmodel.AOTextureMap-1,10);
    case('SSS'); P = TextureSmoothing(Gmodel.Element,Gmodel.SSSTextureMap,10);
    end
    
    set(Gmodel.FigHandle,'FaceVertexCData',P,'facecolor','interp');
    colormap(Gmodel.Colormap);
    drawnow;
end

%--------------------------------------------------------------------- bake
function Gmodel = bake(Gmodel)
    
    Gmodel = BakeCubemap(Gmodel,Gmodel.Texture);
    
%     if Gmodel.SSS, Gmodel.AOInvert = true; end
%     if Gmodel.AO || Gmodel.SSS
%     Gmodel = BakeAmbientOcclusion(Gmodel);
%     end
    if Gmodel.SSS
        Gmodel.SSSRenderComplete = false;
        Gmodel.AOInvert = true;
        Gmodel = BakeAmbientOcclusion(Gmodel);
    end
    
    if Gmodel.AO
        Gmodel.AORenderComplete = false;
        Gmodel.AOInvert = false;
        Gmodel = BakeAmbientOcclusion(Gmodel);
    end
        

end

%---------------------------------------------------------------- export
function export(Gmodel,filename,type)
if nargin < 3, type = 'stl'; end
if nargin < 2
    filename = string(['gmodel','_', char(datetime(now,...
                'ConvertFrom','datenum')),'.stl']);
    filename = erase(filename,[":"," "]);
    type = 'stl'; 
end

fv = struct;
f = Gmodel.Element;
v = Gmodel.Node;

[v, ~, indexn] =  unique(v, 'rows');
f = indexn(f);

if strcmp(type,'stl')
fv.vertices = v;
fv.faces = f;
stlwriter(char(filename),fv);
elseif strcmp(type,'obj')
fv.vertices = v;
fv.faces = f;
fv.objects(1).type='f';
fv.objects(1).data.vertices=fv.faces;
objwriter(fv,'test.obj');

end


end

%---------------------------------------------------------------- slice
function Gmodel = slice(Gmodel,dim,s)
    
if ~isempty(Gmodel.Slice), delete(Gmodel.Slice); end

v = Gmodel.Node;
I = Gmodel.Center(:,3) > s;
J = -(v(:,3) - s) < 0.01*(Gmodel.BdBox(6) - Gmodel.BdBox(5));

N = 150;

lambda = 0.1*(Gmodel.BdBox(6) - Gmodel.BdBox(5));
v = Gmodel.Node(J,:);
B = boxhull(v(:,1:2),lambda);
x = linspace(B(1),B(2),N);
y = linspace(B(3),B(4),N);

[X,Y] = meshgrid(x,y);

fv = struct; 
fv.faces = Gmodel.Element; 
fv.vertices = Gmodel.Node;
INC = double(inpolyhedron(fv,[X(:),Y(:),X(:)*0 + s]));
fprintf('\n');

A = GaussianFilter(reshape(INC,[N,N]),1);
A(A < 0.5) = nan;

hold on;
colormap(autumn(6));

Gmodel.Slice = surf(X,Y,zeros(N,N) + s,'cData',A,'linestyle','-',...
    'AlphaData',single(~isnan(A)),'facecolor','interp','Edgecolor','interp');

caxis([-15 10]);

Gmodel.Node(J,3) = s;
updateNode(Gmodel);
Gmodel.Element(I,:) = nan;

updateElements(Gmodel);
Gmodel.Element = Gmodel.Element0;
Gmodel.Node = Gmodel.Node0;

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
       
    fprintf(['* Loaded mesh = ']);
    cprintf('hyper', [msh{1}, '\n']);   
       
    elseif isa(msh{1},'function_handle')
       if length(msh{2}) ~= 6
           cout('err','* specify a correct bounding box of size 6 x 1');
       end
       
       Gmodel.SDF = msh{1};
       Gmodel.BdBox = msh{2};
       [X,Y,Z,P] = MeshGridding((1+1e-6)*Gmodel.BdBox+1e-6,Gmodel.Quality);
       P = single(P);
       D = Gmodel.SDF(single(P));
       D = reshape(D(:,end),size(X));
       
       [f,v,~] = MarchingCubes(single(X),single(Y),...
          single(Z),single(D),1e-6);

    elseif length(msh) == 2
       v = msh{1};
       f = msh{2};       
    elseif isa(msh{1},'Fem')
       E = msh{1}.Mesh.get('ElemMat');
       v = msh{1}.Node;
       f_ = num2cell(E(:,1:end-1),2);
       V = cellfun(@(x) quad2tri(x),f_,'UniformOutput',false);
       f = vertcat(V{:});
       
    elseif  isa(msh{1},'Mesh')
       v = msh{1}.Node;
       f = msh{1}.Element;
    elseif  isa(msh{1},'Mmesh')
       myGhostFigure = figure("Visible",false);
       [x,y,z] = tubeplot(msh{1}.Node.',msh{1}.get('WireThickness'));
       h = mesh(x,y,z);
       fv = surf2patch(h);
       v = fv.vertices;
       f = fv.faces;
    end
    
    [vn,fn] = TriangleNormal(v,f);
    
    Gmodel.Node = v;
    Gmodel.Node0 = v;
    Gmodel.Element = f;
    Gmodel.Element0 = f;
    Gmodel.VNormal = -vn;
    Gmodel.Normal = fn;  
    Gmodel.RMatrix = eye(4);
    Gmodel.BdBox = boxhull(Gmodel.Node); 
    
    fprintf(['* Vertices  = ', num2str(length(v)/1e3,4), 'k \n']); 
    pause(0.01);
    fprintf(['* Polycount = ', num2str(length(f)/1e3,4), 'k \n']);
    
    fcell = num2cell(Gmodel.Element,2);
    [~,Gmodel.V2F,Gmodel.F2V] = ElementAdjecency(fcell);
  
    Gmodel.TextureMap = Normal2RGB(Gmodel.Normal);
    
    f = num2cell(Gmodel.Element,2);
    VCell = cellfun(@(E) Gmodel.Node(E,:),f,'UniformOutput',false);
    c = cellfun(@(V) mean(V,1),VCell,'UniformOutput',false);
    Gmodel.Center = vertcat(c{:});
end

%------------------------------------------------------------ bake cubemaps
function Gmodel = BakeCubemap(Gmodel,Cubemap)
      
    Ny = size(Cubemap,1); 
    Nx = size(Cubemap,2);
    
    if ishandle(101)
        Phi = view;
    else, 
        Phi = eye(4); 
    end
    
    Phi = transpose(Phi(1:3,1:3));
    
    if ~isempty(Gmodel.View)
       Phi = Phi*Gmodel.View;
    end
    
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

%-------------------------------------------------------- bake texture maps
function Gmodel = BakeAmbientOcclusion(Gmodel)
Bias = Gmodel.AOBias;
Res = round(2.^Gmodel.AOBit);
BB = BoundingBox(Gmodel.Node); 
if Gmodel.AORenderComplete == false
R = Gmodel.AORadius*mean([abs(BB(2)-BB(1)),abs(BB(4)-BB(3)),...
    abs(BB(6)-BB(5))]);
else
R = Gmodel.SSSRadius*mean([abs(BB(2)-BB(1)),abs(BB(4)-BB(3)),...
    abs(BB(6)-BB(5))]);    
end

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
if Gmodel.AORenderComplete == false
    fprintf('* Baking ambient ccclusion... ');
else
    fprintf('* Baking sub-surface scattering... ');
end
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
AOmap = TextureSmoothing(f,AOmap,5);
end

if Gmodel.AORenderComplete == false
    Gmodel.AOTextureMap = AOmap;
    Gmodel.AORenderComplete = true;
end

if Gmodel.SSSRenderComplete == false
    Gmodel.SSSTextureMap = AOmap;
    Gmodel.SSSRenderComplete = true;
end

end

%----------------------------------------------------- generate groundplane
function Groundplane(Gmodel,gnd)
if nargin < 2
    tmp = Gmodel.BdBox; a = 0.1;
    tmp(1) = Gmodel.BdBox(1)-a*( Gmodel.BdBox(2) - Gmodel.BdBox(1)); 
    tmp(2) = Gmodel.BdBox(2)+a*( Gmodel.BdBox(2) - Gmodel.BdBox(1)); 
    tmp(3) = Gmodel.BdBox(3)-a*( Gmodel.BdBox(4) - Gmodel.BdBox(3)); 
    tmp(4) = Gmodel.BdBox(4)+a*( Gmodel.BdBox(4) - Gmodel.BdBox(3)); 
    tmp(5) = Gmodel.BdBox(5)-0*( Gmodel.BdBox(6) - Gmodel.BdBox(5))-1e-3; 
    tmp(6) = Gmodel.BdBox(6)+a*( Gmodel.BdBox(6) - Gmodel.BdBox(5)); 
else
    tmp = gnd;
end

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

%---------------------------------------------------- generate bounding box
function BoundingBox(Gmodel)

tmp = boxhull(Gmodel.Node);

v = [tmp(1),tmp(3), tmp(5);  
     tmp(2),tmp(3), tmp(5);
     tmp(2),tmp(4), tmp(5);
     tmp(1),tmp(4), tmp(5);
     tmp(1),tmp(3), tmp(6);  
     tmp(2),tmp(3), tmp(6);
     tmp(2),tmp(4), tmp(6);
     tmp(1),tmp(4), tmp(6)];
 
f = [1,2,3,4;
     5,6,7,8;
     1,5,nan,nan;
     2,6,nan,nan;
     3,7,nan,nan;
     4,8,nan,nan];

patch('Faces',f,'Vertices',v,...
    'Linewidth',1.5,'linestyle','-','FaceColor','none',...
    'EdgeColor',col(2));

end
    
end
end

