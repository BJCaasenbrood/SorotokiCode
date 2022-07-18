function add2path
%ADD2PATH loads the current folder (and all subfolders within) to MATLAB's
%   search path. 
%
%   Last edit: July 15, 2022.

addpath(genpath(cd));
end

