function x = pingserver
%PINGSERVER True for active internet connection
%
%   See also: REPO, CDSORO.
%
%   For more information: 
%       https://github.com/BJCaasenbrood/SorotokiCode

%   Copyright 2018-2023 B.J.Caasenbrood 
[~,b] = dos('ping -n 1 www.google.com');
n=strfind(b,'Lost');
n1=b(n+7);
if(n1=='0'), x=true; else, x=false; end
end
