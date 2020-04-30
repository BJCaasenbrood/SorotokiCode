function varargout = cdsoro
%CDSORO Returns installation directory of the soft robotics toolkit
%
%   See also REPO.

%   Copyright 2018-2023 B.J.Caasenbrood 
%path = erase(cd,'\src\__base__');    

fid = fopen('startup.m','rt');
C = textscan(fid, '%s', 1, 'CommentStyle', '!'); 
fclose(fid);
C = C{1};
path = erase(C,'%!INSTALDIR:');
path = path{1};

if nargout > 0, varargout{1} = path;
else, cd(path); end

end

