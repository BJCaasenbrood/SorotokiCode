function y = sigmoid(x,k)
if nargin < 2, k = 0.25; 
else, k =  min(max(k,-1+eps),1-eps); end

x = min(max(x,-1),1);

y = (x - k*x)./(k-2*k*abs(x)+1);

end

