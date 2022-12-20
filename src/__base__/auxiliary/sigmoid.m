%% SIGMOID is an smooth activation mapping described by f: Rn -> [-1,1]. 
% The slope of the activation function can be changed by y = sigmoid(x,k) 
% where k is is a positive constant (default k = 0.5). 

function y = sigmoid(x,k)
if nargin < 2
    k = -0.0; 
else
    k =  min(max(k,-1+eps),1-eps); 
end

x = min(max(x,-1),1);            % ensure x is bounded
y = (x - k*x)./(k-2*k*abs(x)+1); % sigmoid function
end

