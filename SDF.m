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
end
end

