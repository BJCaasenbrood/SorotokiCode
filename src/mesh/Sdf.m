classdef Sdf
     
    properties
        sdf;
        BdBox;
        cmap = turbo;
    end
    
    methods        
%         dV;
%         x,y;
%         X,Y;
        
        function obj = Sdf(fnc,varargin)
            obj.sdf = @(x) [fnc(x),fnc(x)];
            obj.cmap = turbo;
            for ii = 1:2:length(varargin)
                obj.(varargin{ii}) = varargin{ii+1};
            end
        end

        function r = plus(obj1,obj2)
            fnc = @(x) dUnion(obj1.sdf(x),obj2.sdf(x));
            r = Sdf(fnc);
            B = [box2node(obj1.BdBox); box2node(obj2.BdBox)];
            r.BdBox = boxhull(B);
        end
        
        function r = minus(obj1,obj2)
            fnc = @(x) dDiff(obj1.sdf(x),obj2.sdf(x));
            r = Sdf(fnc);
            r.BdBox = obj1.BdBox;
        end
        
        function r = mrdivide(obj1,obj2)
            fnc = @(x) dIntersect(obj1.sdf(x),obj2.sdf(x));
            r = Sdf(fnc);
            B = [box2node(obj1.BdBox); box2node(obj2.BdBox)];
            r.BdBox = boxhull(B);
        end
        
        function r = transpose(obj1)
            fnc = @(x) obj1.sdf([x(:,2),x(:,1)]);
            r = Sdf(fnc);
            B = obj1.BdBox;
            r.BdBox = [B(3), B(4), B(1), B(2)];
        end
    end
        
methods (Access = public)

function d = eval(Sdf,x)
    d = Sdf.sdf(x);
end

function show(Sdf,Quality)
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
        
        figure(101);
        surf(X,Y,reshape(D,[Quality Quality]),'linestyle','none');
        axis equal; hold on;
        contour3(X,Y,reshape(D,[Quality Quality]),[0 0],'linewidth',...
            1.5,'Color','w');
        
        colormap(Sdf.cmap);
        view(0,90);
    else
        z = linspace(Sdf.BdBox(5),Sdf.BdBox(6),Quality);
        [X,Y,Z] = meshgrid(x,y,z);

        D = Sdf.eval([X(:),Y(:),Z(:)]);
        D = D(:,end);
        D(D>0) = NaN;
        figure(101);
        C = reshape(D,[Quality Quality,Quality]);
        scatter3(X(:),Y(:),Z(:),85,C(:),'Marker','.');
        axis equal; hold on;

        colormap(Sdf.cmap);
    end
end

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
        
        figure(101);
        contour3(X,Y,reshape(D,[Quality Quality]),[0 0],'linewidth',...
            2,'Color','k'); hold on;
        
        D(D<=0) = 0;
        D(D>1e-6)  = NaN;
        
        surf(X,Y,reshape(D,[Quality Quality]),'linestyle','none');
        axis equal; 
        
        
        colormap(Sdf.cmap);
        view(0,90);
        axis off;
        caxis([-1 1]);
    end
end

end
    
methods (Access = private)

end
end

