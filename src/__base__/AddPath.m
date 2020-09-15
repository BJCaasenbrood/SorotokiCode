function AddPath(str)
%ADDPATH Summary of this function goes here
%   Detailed explanation goes here
if ismac
    % Code to run on Mac platform
elseif isunix
    % Code to run on Linux platform
    str = strrep(str,'\','/');
    addpath(str);
elseif ispc
    % Code to run on Windows platform
    addpath(str);
else

end
end

