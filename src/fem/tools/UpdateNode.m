function [Node,Node0] = UpdateNode(Fem,U)
Node0 = Fem.get('Node0'); 
Ntmp = Node0;

u = FieldOperator(Fem,U);

Ntmp(:,1) = Node0(:,1) + u(:,1);
Ntmp(:,2) = Node0(:,2) + u(:,2);
if Fem.Dim == 3, Ntmp(:,3) = Node0(:,3) + u(:,3); end

Node = Ntmp;
end