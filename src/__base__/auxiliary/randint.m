function y = randint(d1,d2,n)
%RANDINT(X1,X2,n) generates random integer vector between X1 and X2 of
%length n
if nargin == 2
   n = 1;
end 
r = rand(n,1);
y = ceil(d1 + (d2-d1)*r(:)');
end

