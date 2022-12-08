function y = stepspace(x, stepsize, numstep)
%LINSPACE Linearly spaced vector.
%   LINSPACE(X1, X2) generates a row vector of 100 linearly
%   equally spaced points between X1 and X2.
%
%   LINSPACE(X1, X2, N) generates N points between X1 and X2.
%   For N = 1, LINSPACE returns X2.
%
%   Class support for inputs X1,X2:
%      float: double, single
%
%   See also LOGSPACE, COLON.

%   Copyright 1984-2016 The MathWorks, Inc.

if nargin == 2
    numstep = 100;
else
    numstep = floor(double(numstep));
end
if ~isscalar(x) || ~isscalar(stepsize) || ~isscalar(numstep)
    error(message('MATLAB:linspace:scalarInputs'));
end

n1 = numstep-1;
%c = (d2 - d1).*(n1-1); %check intermediate value for appropriate treatment
% if isinf(c)
%     if isinf(d2 - d1) %opposite signs overflow
%         y = d1 + (d2./n1).*(0:n1) - (d1./n1).*(0:n1);
%     else
%         y = d1 + (0:n1).*((d2 - d1)./n1);
%     end
% else
y = x + (0:n1).*stepsize;
%end

% if ~isempty(y)
%     if st == 0
%         y(:) = d1;
%     else
%         y(1) = d1;
%         y(end) = d2;
%     end
% end
