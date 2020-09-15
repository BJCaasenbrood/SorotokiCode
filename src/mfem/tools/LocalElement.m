%--------------------------------------------------------------------------
function [dC,dR,Cm] = LocalElement(Node,Element)
SubElem = num2cell(Element,2);
dC = cellfun(@(E) diff(Node(E,:)),SubElem,'UniformOutput',false);
Cm = cellfun(@(E) mean(Node(E,:)),SubElem,'UniformOutput',false);
dR = cellfun(@(E) norm(E,2),dC,'UniformOutput',false);
%find the relative segment components

dC = vertcat(dC{:}); 
Cm = vertcat(Cm{:}); 
dR = vertcat(dR{:}); 
end