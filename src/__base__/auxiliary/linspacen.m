function C = linspacen(A,B,n)

N = size(A);

if ~all(size(A)==size(B)), error('A and B should be the same size'); end

if n==1
    C=B;
else
    logicReshape=~isvector(A);
   
    A = A(:);
    B = B(:);
    
    C=repmat(A,[1,n])+((B-A)./(n-1))*(0:1:n-1);
    
    %Force start and end to match (also avoids numerical precission issues) 
    C(:,1)=A; %Overide start
    C(:,end)=B; %Overide end
    if logicReshape
        C=reshape(C,[N n]);
    end
end