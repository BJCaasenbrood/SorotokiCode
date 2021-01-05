classdef SDF
    %SDF Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        sdf;
        BdBox;
        cmap = viridis;
    end
    
    methods
        function obj = SDF(fnc)
            obj.sdf = fnc;
        end

        function r = plus(obj1,obj2)
            fnc = @(x) dUnion(obj1.sdf(x),obj2.sdf(x));
            r = SDF(fnc);
            B = [box2node(obj1.BdBox); box2node(obj2.BdBox)];
            r.BdBox = boxhull(B);
        end
        
        function r = minus(obj1,obj2)
            fnc = @(x) dDiff(obj1.sdf(x),obj2.sdf(x));
            r = SDF(fnc);
            r.BdBox = obj1.BdBox;
        end
        
        function r = mrdivide(obj1,obj2)
            fnc = @(x) dIntersect(obj1.sdf(x),obj2.sdf(x));
            r = SDF(fnc);
            B = [box2node(obj1.BdBox); box2node(obj2.BdBox)];
            r.BdBox = boxhull(B);
        end
        
        function r = transpose(obj1)
            fnc = @(x) obj1.sdf([x(:,2),x(:,1)]);
            r = SDF(fnc);
            B = obj1.BdBox;
            r.BdBox = [B(3), B(4), B(1), B(2)];
        end
    end
        
methods (Access = public)

function d = eval(SDF,x)
    d = SDF.sdf(x);
end

function show(SDF,Quality)
    if nargin < 2,
        Quality = 150;
    end
    
    x = linspace(SDF.BdBox(1),SDF.BdBox(2),Quality);
    y = linspace(SDF.BdBox(3),SDF.BdBox(4),Quality);
    [X,Y] = meshgrid(x,y);
    
    D = SDF.eval([X(:),Y(:)]);
%     D = D(:,end);
%     D(D>0) = NaN;
%     
    surf(X,Y,reshape(D(:,end),[Quality Quality]),'linestyle','none');
    axis equal;
    axis equal; hold on;
    contour3(X,Y,reshape(D(:,end),[Quality Quality]),[0 0],'linewidth',...
        1.5,'Color','w');
    colormap(SDF.cmap);
end
end
    
methods (Access = private)

end
end

