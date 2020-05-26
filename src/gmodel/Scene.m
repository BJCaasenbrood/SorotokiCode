classdef Scene

    properties (Access = public)
    Name;
    Object;
    RefNode;
    Build;
    Time;
    State;
    BdBox;
    end
    
    properties (Access = private)
    Figure;
    
    end
    
%--------------------------------------------------------------------------
methods  
    
%------------------------------------------------------------ Process Class
function obj = Scene(varargin) 
    obj.BdBox = [0,0,0,0,0,0];
    obj.Time = 0;
    obj.State = 0;
    for ii = 1:2:length(varargin)
        obj.(varargin{ii}) = varargin{ii+1};
    end
end

%---------------------------------------------------------------------- get     
function varargout = get(Scene,varargin)
    if nargin > 1
        varargout{nargin-1,1} = [];
        for ii = 1:length(varargin)
            varargout{ii,1} = Scene.(varargin{ii});
        end
    else
        varargout = Scene.(varargin);
    end
end
        
%---------------------------------------------------------------------- set
function Scene = set(Scene,varargin)
    for ii = 1:2:length(varargin)
        Scene.(varargin{ii}) = varargin{ii+1};
    end
end

%---------------------------------------------------------------------- set
function Scene = add(Scene,varargin)
    for ii = 1:length(varargin)
        tmp = varargin{ii};
        Scene.Object{end+1,1} = tmp.bake();
        B0 = Scene.BdBox;
        B1 = varargin{ii}.get('BdBox');
        Scene.BdBox = [min(B1(1),B0(1)),max(B1(2),B0(2)),...
            min(B1(3),B0(3)),max(B1(4),B0(4))];
    end
end

%---------------------------------------------------------------------- set
function Scene = render(Scene,varargin)
    
    H = figure(101);
    
    h = rotate3d;
    h.Enable = 'on';
    hold on;
    
    set(gcf,'color',[255, 255, 255]/255); material dull;
    axis equal; axis(Scene.BdBox); axis off;
    daspect([1,1,1]);
    
    ax = gca; ax.Clipping = 'off';
    view(30,30);
    
    for ii = 1:length(Scene.Object)
        gm = Scene.Object{ii};
        
        if strcmp(gm.get('Shading'),'Face'), shading = 'flat';
        else, shading = 'interp';
        end
        
        hp = patch('Vertices',gm.Node,'Faces',gm.Element,'linestyle',...
            'none','FaceVertexCData',gm.get('TextureMap'),'FaceColor',shading);
        
        gm.set('FigHandle',hp);
        gm.update();
        Scene.Object{ii} = gm;
    end
    
    h.ActionPostCallback  = @myprecallback;
    Scene.Figure = H;
    Scene.Figure.UserData = Scene;
    
    function myprecallback(src,evnt)
        update(src.UserData);
    end

end

%--------------------------------------------------------------------------
function Scene = update(Scene)
    
    if ~isempty(Scene.Build)
        for ii = 1:length(Scene.Object)
            Scene.Object{ii} = Scene.Object{ii}.resetNode();
        end
        Scene.Object = Scene.Build(Scene.Object,Scene.State,Scene.Time);
    end
    
    for ii = 1:length(Scene.Object)
        var = Scene.Object{ii}.update();
        set(var{1},'FaceVertexCData',var{2},'Facecolor',var{3},...
            'Vertices',var{4});
    end
    
    axis(Scene.BdBox);
    drawnow;
end

end
end