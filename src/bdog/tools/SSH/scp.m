function ssh2_struct = scp(ssh2_struct)
% SCP   Main SCP interface for SSH2. The input struct provides all 
%       the parameters for the SSH2 function to perform a scp copy.
%  
%   SCP(SSH2_CONN) enables a ssh2 connection for scp and runs
%   it through the ssh2 engine where a previous configured scp command
%   will be executed.
% 
%   SCP returns the SSH2_CONN for future use.
%
%see also scp_get, scp_put, scp_simple_get, scp_simple_put
%
% (c)2011 Boston University - ECE
%    David Scott Freedman (dfreedma@bu.edu)
%    Version 2.0

if nargin == 1
    ssh2_struct.scp = 1;
    ssh2_struct = ssh2(ssh2_struct);
else
    ssh2_struct = [];
    help scp
end