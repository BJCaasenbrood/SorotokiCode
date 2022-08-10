name = 'color_green';
hex = '#2e8540';

HSV = rgb2hsv(hex2rgb(hex));

set = {'darkest','dark','light','lighter','lightest'};
SV  = [0.8235,0.6000;
       0.8434,0.7765;
       0.7309,0.9765;
       0.5412,1.0000;
       0.1804,1.0000];

for ii = 1:numel(set)
    NAME = [name,'_',set{ii}];
    FILE = [NAME,'.m'];
    color = [HSV(1),SV(ii,:)];
    HEX = rgb2hex(hsv2rgb(color));
    writefile(FILE,NAME,HEX);
end

% FID = fopen(file,'W');
% fwrite(FID,['function [y, hex] = ',name,newline]); 
% fwrite(FID,strcat(" hex = ",hex,";"));
% fwrite(FID,newline);
% fwrite(FID,[' y = hex2rgb(hex);',newline]); 
% fwrite(FID,['end']);
% 
% fclose(FID);
% color_secondary_lightest
% hex = '#f9e0de';
% 
% end

function writefile(file,name,hex)
FID = fopen(file,'W');
fwrite(FID,['function [y, hex] = ',name,newline]);
fwrite(FID,strcat(" hex = ","'",hex,"'",";"));
fwrite(FID,newline);
fwrite(FID,[' y = hex2rgb(hex);',newline]);
fwrite(FID,['end']);
fclose(FID);
end