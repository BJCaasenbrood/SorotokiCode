function [Node,Node0] = UpdateNode(Fem,U)
Node0 = Fem.get('Node0'); 
Ntmp = Node0;

[ux,uy,~] = FieldOperator(Fem,U);

Ntmp(:,1) = Node0(:,1) + ux(:);
Ntmp(:,2) = Node0(:,2) + uy(:);

Node = Ntmp;
end

function [ux,uy,un] = FieldOperator(Fem,U)
ux = zeros(Fem.NNode,1);
uy = zeros(Fem.NNode,1);
un = zeros(Fem.NNode,1);

for node = 1:Fem.NNode
    ux(node) = U(2*node - 1,1);
    uy(node) = U(2*node,1);
    un(node) = sqrt(ux(node)^2 + uy(node)^2);
end

end