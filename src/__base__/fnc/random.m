function y = random(d1,d2,n)
%RANDOM(X1,X2,n) generates random integer vector between X1 and X2 of
%length n. 
if nargin == 2
   n = 1;
end 
r = rand(n,1);
y = d1 + (d2-d1)*r(:)';
end

