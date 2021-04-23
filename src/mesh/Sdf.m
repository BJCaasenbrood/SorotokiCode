classdef Sdf
    %SDF Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        sdf;
        BdBox;
        cmap = viridis;
    end
    
    methods
        function obj = Sdf(fnc)
            obj.sdf = @(x) [fnc(x),fnc(x)];
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
    if nargin < 2,
        Quality = 250;
    end
    
    x = linspace(Sdf.BdBox(1),Sdf.BdBox(2),Quality);
    y = linspace(Sdf.BdBox(3),Sdf.BdBox(4),Quality);
    [X,Y] = meshgrid(x,y);
    
    D = Sdf.eval([X(:),Y(:)]);
    D = D(:,end);
    D(D>0) = NaN;
    figure(101);
    surf(X,Y,reshape(D,[Quality Quality]),'linestyle','none');
    axis equal;
    axis equal; hold on;
    contour3(X,Y,reshape(D,[Quality Quality]),[0 0],'linewidth',...
        1.5,'Color','w');
    
    colormap(Sdf.cmap);
end
end
    
methods (Access = private)

end
end

