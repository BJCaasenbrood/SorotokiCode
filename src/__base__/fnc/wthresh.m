function x = wthresh(x,sorh,t)
%MATLAB Code Generation Library Function

%   Copyright 1995-2016 The MathWorks, Inc.
%#codegen

coder.internal.prefer_const(sorh);
SOFT = (sorh == 's');
HARD = (sorh == 'h');
coder.internal.assert(SOFT || HARD, ...
    'Wavelet:FunctionArgVal:Invalid_ArgVal');
if SOFT
    if isscalar(t)
        for k = 1:numel(x)
            x(k) = wthreshSoft(x(k),t(1));
        end
    else
        for k = 1:numel(x)
            x(k) = wthreshSoft(x(k),t(k));
        end
    end
else
    if isscalar(t)
        for k = 1:numel(x)
            x(k) = wthreshHard(x(k),t(1));
        end
    else
        for k = 1:numel(x)
            x(k) = wthreshHard(x(k),t(k));
        end
    end
end

%--------------------------------------------------------------------------

function x = wthreshSoft(x,t)
coder.inline('always');
tmp = abs(x) - t;
tmp = (tmp + abs(tmp))/2;
x = sign(x)*tmp;

%--------------------------------------------------------------------------

function x = wthreshHard(x,t)
coder.inline('always');
x = x*(abs(x) > t);

%--------------------------------------------------------------------------
