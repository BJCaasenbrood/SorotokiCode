function y = random(d1,d2,n)

if nargin == 2
   n = 1;
end 

r = rand(n,1);
y = d1 + (d2-d1)*r(:)';

end

