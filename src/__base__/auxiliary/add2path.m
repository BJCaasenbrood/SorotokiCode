function add2path(varargin)
%ADD2PATH loads the current folder (and all subfolders within) to MATLAB's
%   search path. 
%
%   Last edit: July 15, 2022.
if isempty(varargin)
    addpath(genpath(cd));
else
     addpath(genpath(varargin{1}));
end
end

