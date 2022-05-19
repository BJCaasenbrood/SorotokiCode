function y = stepspace(d1, st, n)
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
    n = 100;
else
    n = floor(double(n));
end
if ~isscalar(d1) || ~isscalar(st) || ~isscalar(n)
    error(message('MATLAB:linspace:scalarInputs'));
end

n1 = n-1;
%c = (d2 - d1).*(n1-1); %check intermediate value for appropriate treatment
% if isinf(c)
%     if isinf(d2 - d1) %opposite signs overflow
%         y = d1 + (d2./n1).*(0:n1) - (d1./n1).*(0:n1);
%     else
%         y = d1 + (0:n1).*((d2 - d1)./n1);
%     end
% else
    y = d1 + (0:n1).*st;
%end

% if ~isempty(y)
%     if st == 0
%         y(:) = d1;
%     else
%         y(1) = d1;
%         y(end) = d2;
%     end
% end
