function [Pc,A,Nv,Ev] = computeCentroid(Mesh)
    
N  = Mesh.NElem;
v  = Mesh.Node;
f  = Mesh.Element;
Pc = zeros(N,2); 
A  = zeros(N,1);
Nv = cell(N,1);
Ev = cell(N,1);

if nargout < 2, cmpNorm = true;
else, cmpNorm = false;
end

for el = 1:N

  vx = v(f{el},1); 
  vy = v(f{el},2); 
  nv = length(f{el}); 
  
  vxS = vx([2:nv 1]); 
  vyS = vy([2:nv 1]); 

  if cmpNorm
      n = ([vx,vy] - [vxS,vyS]);
      eG = zeros(length(vx),1);
      
      for ii = 1:length(n)
          nspace = transpose(null(n(ii,:)));
          if size(nspace,1) == 2, n(ii,:) = nspace(1,:);
          else, n(ii,:) = nspace;
          end
          a = [vxS(ii)-vx(ii),vyS(ii)-vy(ii)];
          eG(ii,1) = sqrt(a(1)^2 + a(2)^2);
          a = a/eG(ii,1);
          b = [n(ii,1),n(ii,2)];
          d = b(1)*a(2) - b(2)*a(1);
          n(ii,:) = sign(d)*n(ii,:);
      end
      
      nS = n([nv 1:nv-1],:); n = nS+n;
      eS = eG([nv 1:nv-1],:);
      Nv{el} = (n./sqrt(n(:,1).^2 + n(:,2).^2));
      Ev{el} = 0.5*(eS + eG);

  end
  
  tmp = vx.*vyS - vy.*vxS;
  A(el) = 0.5*sum(tmp);
  Pc(el,:) = 1/(6*A(el,1))*[sum((vx+vxS).*tmp),...
                            sum((vy+vyS).*tmp)];
end

end