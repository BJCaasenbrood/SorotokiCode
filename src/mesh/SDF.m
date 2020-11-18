classdef SDF
    %SDF Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        sdf
    end
    
    methods
        function obj = SDF(fnc)
            obj.sdf = fnc;
        end

        function r = plus(obj1,obj2)
            fnc = @(x) dUnion(obj1.sdf(x),obj2.sdf(x));
            r = SDF(fnc);
        end
        
        function r = minus(obj1,obj2)
            fnc = @(x) dDiff(obj1.sdf(x),obj2.sdf(x));
            r = SDF(fnc);
        end
        
        function r = mrdivide(obj1,obj2)
            fnc = @(x) dIntersect(obj1.sdf(x),obj2.sdf(x));
            r = SDF(fnc);
        end
    end
        
methods (Access = public)

function d = eval(SDF,x)
    d = SDF.sdf(x);
end
end
    
methods (Access = private)
function d = dUnion(d1,d2) 
    d=[d1(:,1:(end-1)),d2(:,1:(end-1))];
    d=[d,min(d1(:,end),d2(:,end))];
end

function d = dDiff(d1,d2) 
d=[d1(:,1:(end-1)),d2(:,1:(end-1))];
d=[d,max(d1(:,end),-d2(:,end))];
end

function d = dIntersect(d1,d2) % max(d1,d2)
d=[d1(:,1:(end-1)),d2(:,1:(end-1))];
d=[d,max(d1(:,end),d2(:,end))];
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



end
end

