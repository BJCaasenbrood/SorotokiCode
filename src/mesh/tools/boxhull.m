function BdBox = boxhull(Node,eps)
if nargin < 2, eps = 1e-6; end

if size(Node,2) == 2
    BdBox = zeros(4,1);
    BdBox(1) = min(Node(:,1))-eps;
    BdBox(2) = max(Node(:,1))+eps;
    BdBox(3) = min(Node(:,2))-eps;
    BdBox(4) = max(Node(:,2))-eps;
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

