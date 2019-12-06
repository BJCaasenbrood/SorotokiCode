function WriteToFile(str)
global FID;
str = strrep(str,'\','\\');
fprintf(FID, ['addpath(genpath(''',str, '''))\n']);
pause(0.01);
cprintf('Keyword',[str,'\n']);
end
