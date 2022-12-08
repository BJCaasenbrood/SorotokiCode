function y = sawtooth(t,varargin)
%SAWTOOTH y = sawtooth(t)
if isempty(varargin)
    window = 1;
else
    window = varargin{1};
end
    y = mod(t,window)/window;
end

