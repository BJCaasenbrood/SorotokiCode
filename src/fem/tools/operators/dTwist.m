function d = dTwist(d1,k) 
dx = d1(:,1); DX = num2cell(dx);
dy = d1(:,2); DY = num2cell(dy);
dz = d1(:,3); DZ = num2cell(dz);

dmax = max(dz); 
dmin = min(dz); 
scale = k*(2*pi)/(dmax-dmin);

dtwist = cellfun(@(E,V,F) Twistoperation(E,V,F,scale),...
    DX,DY,DZ,'UniformOutput',false);

d=vertcat(dtwist{:});

end
%-------------------------------------------------------------------------%

function d = Twistoperation(dx,dy,dz,k)
R = [cos(k*dz),-sin(k*dz),0;...
     sin(k*dz), cos(k*dz),0;...
             0,         0, 1];

d = R*[dx;dy;dz];
d = d(:).';
end