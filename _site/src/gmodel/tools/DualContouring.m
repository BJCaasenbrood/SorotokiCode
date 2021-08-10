function [F,V] = DualContouring(x,y,z,c,n,iso)
F = []; V = [];

n = size(c) - 1;
cc = zeros(n(1),n(2),n(3),'uint16');

vertex_idx = {1:n(1), 1:n(2), 1:n(3); ...
              2:n(1)+1, 1:n(2), 1:n(3); ...
              2:n(1)+1, 2:n(2)+1, 1:n(3); ...
              1:n(1), 2:n(2)+1, 1:n(3); ...
              1:n(1), 1:n(2), 2:n(3)+1; ...
              2:n(1)+1, 1:n(2), 2:n(3)+1; ...
              2:n(1)+1, 2:n(2)+1, 2:n(3)+1; ...
              1:n(1), 2:n(2)+1, 2:n(3)+1 };

for ii=1:8                             % loop thru vertices of all cubes
    idx = c(vertex_idx{ii, :}) > iso;  % which cubes have vtx ii > iso
    cc(idx) = bitset(cc(idx), ii);     % for those cubes, turn bit ii on
end


end

function v = dualContourFitVertex(f,f_normal,x,y,z)

if ~ADPATIVE
   v = [x+DX; y+DX; z+DX];
end

v = zeros(2,2,2);

for dx = [0,DX]
    for dy = [0,DX]
        for dz = [0,DX]
           v() = f([x+dx,y+dy,z+dz]);
        end
    end
end

end