function u = FieldOperator(Fem,U)
mm = Fem.Dim;
u = zeros(Fem.NNode,mm);

for node = 1:Fem.NNode
   u(node,1) = U(mm*node + (1-mm),1);
   u(node,2) = U(mm*node + (2-mm),1);
   if Fem.Dim == 3, u(node,3) = U(mm*node,1); end
end

end