function y = smoothpulse(t,t1,t2,varargin)
if isempty(varargin)
    a = 1;
else
    a = 1/clamp(varargin{1},1e-3,1);
end

y = smoothstep(a*(t-t1)).*(1 - smoothstep(a*(t-t2)));
end

