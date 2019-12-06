%------------------------------------- GENERATE ADJECENCY MATRIX FROM FACES
function [A,F2V,V2F] = ElementAdjecency(face)
n = length(face);    
p = max(cellfun(@max,face));    
MaxNVer = max(cellfun(@numel,face));     
tri = GenerateElementalMatrix(face);
set = 1:n;

countsperrow = SetupM(tri);
idx = 0;
I = zeros(sum(countsperrow),1);
J = zeros(sum(countsperrow),1);
for ii = 1:MaxNVer
   id = ~isnan(tri(:,ii));
   i = set(id)';
   j = tri(id,ii);
   I(idx+1:idx+countsperrow(ii)) = i;
   J(idx+1:idx+countsperrow(ii)) = j;
   idx = idx + countsperrow(ii);
end

M = sparse(I,J,1,n,p);

C = M*M';
C=C-diag(diag(C));

A = sparse(double(C>1));
V2F = (M./sum(M,1))';
F2V = ((M')./sum(M,2)')';
end

%------------------------------------------------ GENERATE ELEMENTAL MATRIX
function cnt = SetupM(tri)
maxver = num2cell(1:size(tri,2));
fnc = @(x) sum(tri(:,x)>=1);
cnt = cellfun(@(x) fnc(x), maxver);
end

%------------------------------------------------ GENERATE ELEMENTAL MATRIX
function A = GenerateElementalMatrix(face)
Element = face(1:length(face))';                 
MaxNVer = max(cellfun(@numel,Element));      
PadWNaN = @(E) [E NaN(1,MaxNVer-numel(E))]; 
ElemMat = cellfun(PadWNaN,Element,'UniformOutput',false);
ElemMat = vertcat(ElemMat{:});       
A = ElemMat;
end
