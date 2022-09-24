classdef Sdf
     
    properties
        sdf;
        BdBox;
        cmap;
        eps = 1e-5;
        
        Node;      
        Element;
        Sample;
        Quality;
    end
    
%--------------------------------------------------------------------------    
    methods        
%---------------------------------------------------- Signed Distance Class     
        function obj = Sdf(fnc,varargin)
            obj.sdf = @(x) [fnc(x),fnc(x)];
            obj.cmap = redgreen;
            obj.Quality = 50;
            for ii = 1:2:length(varargin)
                obj.(varargin{ii}) = varargin{ii+1};
            end
        end
%-------------------------------------------------------------------- union
        function r = plus(obj1,obj2)
            fnc = @(x) dUnion(obj1.sdf(x),obj2.sdf(x));
            r = Sdf(fnc);
            B1 = box2node(obj1.BdBox); 
            B2 = box2node(obj2.BdBox);
            
            if norm(B1) == Inf
                r.BdBox = B2;
            elseif norm(B2) == Inf
                r.BdBox = B1;
            else
                B = [box2node(obj1.BdBox); box2node(obj2.BdBox)];
                r.BdBox = boxhull(B);
            end
        end
%--------------------------------------------------------------- difference        
        function r = minus(obj1,varargin)
            obj2 = varargin{1};
            fnc = @(x) dDiff(obj1.sdf(x),obj2.sdf(x));
            r = Sdf(fnc);
            r.BdBox = obj1.BdBox;
        end
%---------------------------------------------------------------- intersect        
        function r = mrdivide(obj1,obj2)
            fnc = @(x) dIntersect(obj1.sdf(x),obj2.sdf(x));
            r = Sdf(fnc);
            B1 = box2node(obj1.BdBox); 
            B2 = box2node(obj2.BdBox);
            
            if norm(B1) == Inf
                r.BdBox = B2;
            elseif norm(B2) == Inf
                r.BdBox = B1;
            else
                B = [box2node(obj1.BdBox); box2node(obj2.BdBox)];
                r.BdBox = boxhull(B);
            end
        end
%---------------------------------------------------------------- intersect                
        function r = repeat(obj,dX,varargin)
            if isempty(varargin),
                N = 2;
            else
                N = varargin{1};
            end
                
            fnc = @(x) obj.sdf(pRepeat(x,dX));
            r = Sdf(fnc);
            B = obj.BdBox;
            DDX = B(2) - B(1); 
            DDY = B(4) - B(3);
            DDN = norm(dX);
            if dX(2) == 0
                r.BdBox = (obj.BdBox(:) + [0,(N-1)*(DDX),0,0].').';
            elseif dX(1) == 0
                r.BdBox = (obj.BdBox(:) + [0,0,0,(N-1)*(DDY)].').';
            end
        end
%----------------------------------------------------------- rotate X <-> Y              
        function r = transpose(obj1)
            
            B = obj1.BdBox;
            if numel(B) == 4
                fnc = @(x) obj1.sdf([x(:,2),x(:,1)]);
                r = Sdf(fnc);
                r.BdBox = [B(3), B(4), B(1), B(2)];
                r.Node = [obj1.Node(:,2), obj1.Node(:,1)];
                r.Element = obj1.Element;
            else
                fnc = @(x) obj1.sdf([x(:,2),x(:,3),x(:,1)]);
                r = Sdf(fnc);
                r.BdBox = [B(5), B(6), B(3), B(4), B(1), B(2)];
%                 r.Node = [obj1.Node(:,2), obj1.Node(:,3), obj1.Node(:,1)];
%                 r.BdBox = boxhull(r.Node);
%                 r.Element = obj1.Element;
            end
            
        end
%----------------------------------------------------------- rotate X <-> Y          
        function r = translate(obj1,move)
            
           B = obj1.BdBox; 
           if numel(move) == 2
               y = @(x) x - repmat(move(:)',length(x),1);
               fnc = @(x) obj1.sdf(y(x));
               r = Sdf(fnc);
               r.BdBox = [B(1)+move(1), B(2)+move(2),...
                          B(3)+move(1), B(4)+move(2)];
               r.Node = [obj1.Node(:,1)+move(1), obj1.Node(:,2)+move(2)];
               r.Element = obj1.Element;
           else
               y = @(x) x - repmat(move(:)',length(x),1);
               fnc = @(x) obj1.sdf(y(x));
               r = Sdf(fnc);
               r.BdBox = [B(1)+move(1), B(2)+move(1),...
                          B(3)+move(2), B(4)+move(2),...
                          B(5)+move(3), B(6)+move(3)];
%                r.Node = [obj1.Node(:,1)+move(1), obj1.Node(:,2)+move(2),...
%                    obj1.Node(:,3)+move(3)];
%                r.Element = obj1.Element;
           end
           
        end        
%----------------------------------------------------------- rotate X <-> Y              
        function r = rotate(obj1,varargin)
            if isempty(varargin), k = pi/2;
            else, k = varargin{1};
            end
            R = TwoDRotation(k);
            fnc = @(x) obj1.sdf((R*x.').');
            r = Sdf(fnc);
            B = obj1.BdBox;
            BB = [B(1),B(3);B(1),B(4);B(2),B(3);B(2),B(4)];
            
            r.BdBox = boxhull((R*BB.').').';
        end        
%----------------------------------------------------------- rotate X <-> Y              
        function r = smoothunion(obj1,obj2,varargin)
            if isempty(varargin), k = 1;
            else, k = varargin{1};
            end
                
            fnc = @(x) dSmoothUnion(obj1.sdf(x),obj2.sdf(x),k);
            r = Sdf(fnc);
            B1 = box2node(obj1.BdBox); 
            B2 = box2node(obj2.BdBox);
            
            if norm(B1) == Inf
                r.BdBox = B2;
            elseif norm(B2) == Inf
                r.BdBox = B1;
            else
                B = (1+k/5)*[B1;B2];
                r.BdBox = boxhull(B);
            end
        end     
%----------------------------------------------------------- rotate X <-> Y              
        function r = smoothdiff(obj1,obj2,varargin)
            if isempty(varargin), k = 1;
            else, k = varargin{1};
            end
                
            fnc = @(x) dSmoothDiff(obj1.sdf(x),obj2.sdf(x),k);
            r = Sdf(fnc);
            B1 = box2node(obj1.BdBox); 
            B2 = box2node(obj2.BdBox);
            
            if norm(B1) == Inf
                r.BdBox = B2;
            elseif norm(B2) == Inf
                r.BdBox = B1;
            else
                B = [box2node(obj1.BdBox); box2node(obj2.BdBox)];
                r.BdBox = boxhull(B);
            end
        end          
%----------------------------------------------------------- rotate X <-> Y    
        function r = shell(obj,T)
            fnc = @(x) abs(obj.sdf(x)) - T/2;
            r = Sdf(fnc);
            B = obj.BdBox;
            if numel(obj.BdBox) == 6
            r.BdBox = [B(1)-T/2,B(2)+T/2,...
                       B(3)-T/2,B(4)+T/2,...
                       B(5)-T/2,B(6)+T/2];
            else
            r.BdBox = [B(1)-T/2,B(2)+T/2,...
                       B(3)-T/2,B(4)+T/2]
            end
        end
    end
%--------------------------------------------------------------------------       
methods (Access = public)
%--------------------------------------------------------- evalution of SDF
function d = eval(Sdf,x)
d = Sdf.sdf(x);
end
%--------------------------------------------------------- evalution of SDF
function y = intersect(Sdf,x)
d = Sdf.sdf(x);
y = d(:,end)<0;
end
%--------------------------------------------------------- evalution of SDF
function [Nds,X,Y] = sampleSet(Sdf,Quality)
if nargin < 2
    if numel(Sdf.BdBox) < 6, Sdf.Quality = 250;
    else, Sdf.Quality = 50;
    end
end

x = linspace(Sdf.BdBox(1),Sdf.BdBox(2),Sdf.Quality);
y = linspace(Sdf.BdBox(3),Sdf.BdBox(4),Sdf.Quality);
[X,Y] = meshgrid(x,y);
Nds = [X(:), Y(:)];

end
%------------------------------------------------- evalution of tangent SDF
function [T,N,B,Z] = normal(Sdf,x)
d = Sdf.sdf(x);
N = zeros(size(x,1),3);

if size(x,2) == 2
    n1 = (Sdf.sdf(x+repmat([Sdf.eps,0],size(x,1),1))-d)/Sdf.eps;
    n2 = (Sdf.sdf(x+repmat([0,Sdf.eps],size(x,1),1))-d)/Sdf.eps;
    %[Fy] = gradient(d(:,end));   
    %T = [Fy(:,1),Fy(:,2)];
    N(:,1) = n1(:,end);
    N(:,2) = n2(:,end);
else
    n1 = (Sdf.sdf(x+repmat([Sdf.eps,0,0],size(x,1),1))-d)/Sdf.eps;
    n2 = (Sdf.sdf(x+repmat([0,Sdf.eps,0],size(x,1),1))-d)/Sdf.eps;
    n3 = (Sdf.sdf(x+repmat([0,0,Sdf.eps],size(x,1),1))-d)/Sdf.eps;
    N(:,1) = n1(:,end);
    N(:,2) = n2(:,end);
    N(:,3) = n3(:,end);
    %N = [n1(:,end),n2(:,end),n3(:,end)];
end

N = N./sqrt((sum((N.^2),2)));
T = ([0,1,0;-1,0,0;0,0,1]*N.').';
B = cross(T,N);
Z = atan2(N(:,2).*sign(-d(:,end)),N(:,1).*sign(-d(:,end)));
    
end
%--------------------------------------------------------------------- show
function show(Sdf,varargin)
%     
%     if nargin < 2 
%         if numel(Sdf.BdBox) < 6, Quality = 250;
%         else, Quality = 50;
%         end
%     end

    for ii = 1:2:length(varargin)
        Sdf.(varargin{ii}) = varargin{ii+1};
    end
    
    x = linspace(Sdf.BdBox(1),Sdf.BdBox(2),Sdf.Quality);
    y = linspace(Sdf.BdBox(3),Sdf.BdBox(4),Sdf.Quality);
    
    if numel(Sdf.BdBox) < 6
        [X,Y] = meshgrid(x,y);

        D = Sdf.eval([X(:),Y(:)]);
        D = abs(D(:,end)).^(0.5).*sign(D(:,end));
        
        figure(101);
        cplane(X,Y,reshape(D,[Sdf.Quality Sdf.Quality])-1e-6);
        axis equal; hold on;
        
        if isempty(Sdf.Element)
            contour3(X,Y,reshape(D,[Sdf.Quality Sdf.Quality]),[0 0],'linewidth',...
                2.5,'Color','w');
        else
            patch('Vertices',Sdf.Node,'Faces',Sdf.Element,...
                'LineW',2.5,'EdgeColor','w');
        end
        
        cmax = max(abs(D));
        caxis([-cmax,1.25*cmax]);
        axis(Sdf.BdBox);
        colormap(Sdf.cmap);
        view(0,90);
        
    else
        z = linspace(Sdf.BdBox(5),Sdf.BdBox(6),Sdf.Quality);
        [X,Y,Z] = meshgrid(x,y,z);

        D = Sdf.eval([X(:),Y(:),Z(:)]);
        D = D(:,end);
        D(D>0) = NaN;
        figure(101);
        C = reshape(D,[Sdf.Quality Sdf.Quality,Sdf.Quality]);
        scatter3(X(:),Y(:),Z(:),85,C(:),'Marker','.');
        axis equal; hold on;

        colormap(Sdf.cmap);
    end
end

function skeleton(Sdf)
    patch('Vertices',Sdf.Node,'Faces',Sdf.Element,...
                'LineW',2.5,'EdgeColor',col(1),'FaceColor','none');
end

%------------------------------------------------------------- show contour
function showcontour(Sdf,Quality)
    if nargin < 2 
        if numel(Sdf.BdBox) < 6
            Quality = 250;
        else
            Quality = 25;
        end
    end
    
    x = linspace(Sdf.BdBox(1),Sdf.BdBox(2),Quality);
    y = linspace(Sdf.BdBox(3),Sdf.BdBox(4),Quality);
    
    if numel(Sdf.BdBox) < 6
        [X,Y] = meshgrid(x,y);

        D = Sdf.eval([X(:),Y(:)]);
        D = D(:,end);
        
        %figure(101);
        contour3(X,Y,reshape(D,[Quality Quality]),[0 0],'linewidth',...
            2,'Color','k'); hold on;
        
        D(D<=0) = 0;
        D(D>1e-6)  = NaN;
        
        cplane(X,Y,reshape(D,[Quality Quality]));
        %axis equal; 
        
        
        colormap(Sdf.cmap);
        caxis([-1 1]);
    else
        hold on;
        patch('Vertices',Sdf.Node,'Faces',Sdf.Element,...
            'FaceColor','none','LineW',3);
    end
end
end
    
methods (Access = private)

end
end

function [R] = TwoDRotation(x)
R = [cos(x),sin(x);-sin(x),cos(x)];
end
function BdBox = boxhull(Node,eps)
if nargin < 2, eps = 1e-6; end

if size(Node,2) == 2
    BdBox = zeros(4,1);
    BdBox(1) = min(Node(:,1));
    BdBox(2) = max(Node(:,1));
    BdBox(3) = min(Node(:,2));
    BdBox(4) = max(Node(:,2));
else
    BdBox = zeros(6,1);
    BdBox(1) = min(Node(:,1))-eps;
    BdBox(2) = max(Node(:,1))+eps;
    BdBox(3) = min(Node(:,2))-eps;
    BdBox(4) = max(Node(:,2))+eps;
    BdBox(5) = min(Node(:,3))-eps;
    BdBox(6) = max(Node(:,3))+eps;
end

end