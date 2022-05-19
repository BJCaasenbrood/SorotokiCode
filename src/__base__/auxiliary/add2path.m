% function add2path
%-------------------------------------------------------------------------
% Function loads current folder (and all subfolders within) to MATLAB's 
% search path. 
%------------------------------------------------------------------------
function add2path
addpath(genpath(cd));
end

