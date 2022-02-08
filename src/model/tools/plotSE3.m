function h = plotSE3(H,varargin)

if size(H,2) == 7 || size(H,1) == 7
    p = H(5:7);
    R = quat2rot(H(1:4));
else
    p = H(1:3,4);
    R = H(1:3,1:3);
end

if nargin < 2
   a = 1; 
else
   a = varargin{1};
end

for ii = 1:3
    h{ii} = plotarrow(p,a*R(:,ii),ii);
end

end

function h = plotarrow(p,n,id)
hold on; 
h = quiver3(p(1),p(2),p(3),...
    n(1), n(2), n(3),'Color',col(id));
end


