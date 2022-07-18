function x = pingserver
%PINGSERVER returns true for active internet connection. It will simply
%   if www.google.com can be reached. Usage:
%
%   if pingserver
%       disp('Hello internet!');
%   else
%       disp('404')
%   end
%
%   Last edit: Jul 15, 2022.
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
