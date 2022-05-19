function x = pingserver
%PINGSERVER True for active internet connection
%
%   See also: REPO, CDSORO.
%
%   For more information: 
%       https://github.com/BJCaasenbrood/SorotokiCode

%   Copyright 2018-2023 B.J.Caasenbrood 
if ismac
    % Code to run on Mac platform
elseif isunix
    % Code to run on Linux platform
    [~,b] = system('ping -c 1 www.google.com');
    n=strfind(b,'loss');
    n1=b(n-10);
elseif ispc
    [~,b] = system('ping -n 1 www.google.com');
    n=strfind(b,'Lost');
    n1=b(n+7);
else
    disp('Platform not supported')
end

if(n1=='0'), x=true; else, x=false; end

end
