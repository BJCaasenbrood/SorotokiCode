clr;
%% make signed distance func.
sdf = sRectangle(10,10);

%% make mesh
msh = Mesh(sdf,'Quads',[9,9]);
msh = msh.generate();

%% show mesh;
I = [11; 12; 13; 15; 16; 17;
    31; 32; 33; 38; 39; 40; 41; 42; 43; 44; 
    49; 50; 51;
    65; 66; 67; 69; 70; 71];

% I = [14; 23;
%     31; 32; 33; 38; 39; 40; 41; 42; 43; 44; 
%     49; 50; 51;
%     59; 68];

% I = [
%     31; 32; 33; 40; 41; 42; 
%     49; 50; 51;];

% msh = msh.remove(I);
% msh = msh.rebuild();

msh.show();