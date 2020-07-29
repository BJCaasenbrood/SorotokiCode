function ssh2_struct = ssh2(ssh2_struct, arg2, arg3, arg4, arg5)
% SSH2   main control program for interfacing with Ganymed Java SSH-2 Library
%
% SSH2(SSH2_CONN) uses the SSH2_CONN to perform a SSH2 function, e.g.
% issue a command or transfer a file using SCP or SFTP. All options for the
% SSH2 connection are specified in the SSH2_CONN struct. This struct can be
% configued manually, or by using any of the SSH2 support functions
% available. 
%
% The SSH2_CONN is returned, but can be reused. Typically, the commands and
% scp/sftp options are cleared after each function call. If you don't
% desire this behaviour, supply a copy of the SSH2_CONN into the SSH2
% function.
%
% The connection is not closed be default. To close a connection, use
% SSH2_CLOSE, which will setup SSH2_CONN to close and then call this
% function. The exception to this are the SIMPLE functions. They all close
% the connection be default.
%
% For convience, SSH2 also provides a front-end to the SSH2_SIMPLE_COMMAND. 
% [SSH2_STRUCT] = SSH2(HOSTNAME,USERNAME,PASSWORD,COMMAND,[ENABLEPRINTTOSCREEN])
%
%   SSH2 returns the SSH2_CONN for future use.
%
%see also ssh2_config, ssh2_config_publickey, ssh2_command, scp_get, scp_put, sftp,get, sftp_put, ssh2_close
%
% (c)2011 Boston University - ECE
%    David Scott Freedman (dfreedma@bu.edu)
%    Version 2.0

%% BEGIN CODE

input_check = 0; % have not passed the input check
if nargin < 1
    ssh2_struct = [];
else
    if (isstruct(ssh2_struct))
        if (isfield(ssh2_struct,'verified_config') || ...
           isfield(ssh2_struct,'hostname') && ...
           isfield(ssh2_struct,'username') && ...
           isfield(ssh2_struct,'password') ) % min required inputs
                input_check = 1;
                ssh2_struct.verified_config = 1;
                if (nargin == 2)
                    ssh2_struct.command = arg2;
                end
        else

        end
    else
        if nargin >= 4
            input_check = 2; %alternate use of ssh2
        end
    end
end

% do we have enough to proceed?
if (input_check == 0)
    help ssh2.m;
    %fprintf('\n\n!!Incorrect input supplied to ssh2, please try again!!\n');
    if nargout == 0
        clear ssh2_struct
    end
elseif (input_check == 2) %allow ssh2 to be interface to ssh_simpleCommand
    hostname = ssh2_struct;
    username = arg2;
    password = arg3;
    command = arg4;
    if nargin < 5
        if nargout == 0
            enableprint = 1;
        else
            enableprint = 0;
        end
    else
        enableprint = arg4;
    end
    if nargout == 0
        ssh2_simple_command(hostname, username, password, command, enableprint);
    else
        ssh2_struct = ssh2_simple_command(hostname, username, password, command, enableprint);
    end
else
    ssh2_struct = ssh2_main(ssh2_struct);    
end