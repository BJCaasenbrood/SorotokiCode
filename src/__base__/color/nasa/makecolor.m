name = 'color_visited';
hex = "'#4c2c92'";
file = [name,'.m'];

FID = fopen(file,'W');
fwrite(FID,['function [y, hex] = ',name,newline]); 
fwrite(FID,strcat(" hex = ",hex,";"));
fwrite(FID,newline);
fwrite(FID,[' y = hex2rgb(hex);',newline]); 
fwrite(FID,['end']);

fclose(FID);
% color_secondary_lightest
% hex = '#f9e0de';
% 
% end
