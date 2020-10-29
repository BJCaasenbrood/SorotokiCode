classdef Rig < handle

    properties (Access = public)
        Model;
        List;
        IKlist;
        Frame;
        g;
        g0;
    end
    
    properties (Access = private)
        SweepHandle;
        Domain;
        ListDraw;
        AutoScale;
        RigHandle;
    end
    
%--------------------------------------------------------------------------
methods  
%---------------------------------------------------------------- Fem Class
function obj = Rig(mdl,varargin) 
    
    obj.Model = mdl;
    obj.Domain = mdl.get('Sdomain');
    obj.AutoScale = true;
    obj.Frame = 1;
    obj.g0 = [1,0,0,0,0,0,0];
    
    for ii = 1:2:length(varargin)
        obj.(varargin{ii}) = varargin{ii+1};
    end
end

%---------------------------------------------------------------------- get     
function varargout = get(Rig,varargin)
    if nargin > 1
        varargout{nargin-1,1} = [];
        for ii = 1:length(varargin)
            varargout{ii,1} = Rig.(varargin{ii});
        end
    else
        varargout = Rig.(varargin);
    end
end
        
%---------------------------------------------------------------------- set
function Rig = set(Rig,varargin)
    for ii = 1:2:length(varargin)
        Rig.(varargin{ii}) = varargin{ii+1};
    end
end

%---------------------------------------------------------------------- add
function Rig = add(Rig,varargin)
    
    for ii = 1:length(varargin)
        
        gmdl = varargin{ii};
        
        if Rig.AutoScale
           gmdl = gmdl.set('Node0',gmdl.Node*Rig.Domain);
           gmdl = gmdl.set('Node',gmdl.Node*Rig.Domain); 
           gmdl.Texture = grey;
           gmdl = gmdl.bake();
        end

        Rig.List{end+1,1} = gmdl;
        Rig.ListDraw{end+1} = true;
    end
end

%---------------------------------------------------------------------- add
function Rig = reset(Rig)
    for ii = 1:length(Rig.List)
        Rig.List{ii} = Rig.List{ii}.reset();
    end
end

%---------------------------------------------------------------------- rig
function Rig = parent(Rig,Child,ChildNode,ParentNode)
    Rig.IKlist(end+1,:) = [Child,ChildNode,ParentNode];
end

%---------------------------------------------------------------------- rig
function Rig = scale(Rig,Child,Scale)
    gmdl = Rig.List{Child};
    gmdl = Blender(gmdl,'Scale',Scale);
    gmdl = gmdl.set('Node0',gmdl.Node);
    gmdl = gmdl.bake();
    Rig.List{Child} = gmdl;
end

%---------------------------------------------------------------------- rig
function Rig = texture(Rig,Child,TexMap)
    for ii = 1:length(Child)
        Rig.List{Child(ii)}.Texture = TexMap;
    end
end

%---------------------------------------------------------------------- rig
function Rig = hide(Rig,varargin)
    for ii = 1:length(varargin)
        Rig.ListDraw{varargin{ii}} = false;
    end
end

%----------------------------------------------------------- compute ik rig
function Rig = compute(Rig,t)
    
    Rig = Rig.reset();
    [Rig.g,X] = Rig.Model.string(t,200);
    
    s = Rig.Domain;
    N = length(Rig.List);
    Instr = cell(1,1);
    
    for ii = 1:N
        [I,~] = find(Rig.IKlist(:,1) == ii);
        if numel(I) == 1, 
            Instr{ii,1} = 'SE3';

            [~,id] = min(abs(X - s*Rig.IKlist(I,3)));
            Instr{ii,2} = id;

        else, 
            Instr{ii,1} = 'Sweep';
            [a,~] = find(Rig.IKlist(I,2) == 0);
            [b,~] = find(Rig.IKlist(I,2) == 1);
            v = Rig.IKlist(I,3);
            
            id = find(X >= s*v(a) & X <= s*v(b));
            Instr{ii,2} = id;

            if (v(b)-v(a)) < 1.0 && Rig.AutoScale
                gmdl = Rig.List{ii};
                gmdl = Blender(gmdl,'Scale',(v(b)-v(a)));
                %Rig.List{ii} = gmdl.set('Node0',gmdl.Node);
            end
            
            Rig.List{ii} = Blender(Rig.List{ii},'Translate',{'z',v(a)*Rig.Domain});
        end
    end
    
    for ii = 1:size(Instr,1)
        if strcmp(Instr{ii,1},'Sweep')
           id = Instr{ii,2};
           LinkID = knnsearch(X(id).',Rig.List{ii}.Node(:,3));
           Rig.List{ii} = Blender(Rig.List{ii},'Sweep', {LinkID,Rig.g(id,:)});
        elseif strcmp(Instr{ii,1},'SE3')
           Rig.List{ii} = Blender(Rig.List{ii},'SE3', Rig.g(Instr{ii,2},:));
        end
    end
    
    for ii = 1:length(Rig.List)
        Rig.List{ii} = Blender(Rig.List{ii},'SE3', Rig.g0);
    end
       
end

%---------------------------------------------------------------------- rig
function Rig = render(Rig) 
    
   for ii = 1:length(Rig.List)
       if Rig.ListDraw{ii}
        Rig.List{ii}.bake().render();
       end
   end
   
end

%---------------------------------------------------------------------- rig
function Rig = update(Rig) 
    
   for ii = 1:length(Rig.List)
       if Rig.ListDraw{ii}
         Rig.List{ii}.update();
       end
   end
   
end

%---------------------------------------------------------------------- rig
function Rig = showSweep(Rig)
    
    p = Rig.g(:,5:end);
    R = quat2rot(Rig.g0(1:4));
    
    p = (R*p.').' + Rig.g0(5:end);
    
    if ~isempty(Rig.SweepHandle)
       delete(Rig.SweepHandle);
    end
    
    figure(101); hold on;
    Rig.SweepHandle = plot3(p(:,3),p(:,2),p(:,1),...
    'Linewidth',2,'Color',col(2));
    
end

end

methods (Access = private)
    
%----------------------------------------------------- generate groundplane
function Groundplane(Rig,B)
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

