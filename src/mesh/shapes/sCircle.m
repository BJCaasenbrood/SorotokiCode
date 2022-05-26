function sdf = sCircle(xc,yc,r)

if nargin < 2
    r = xc; 
    xc = 0;
    yc = 0;
end
eps = 1e-4*r;
sdf = Sdf(@(P) dCircle(P,xc,yc,r));
sdf.BdBox = [xc-r-eps,xc+r+eps,yc-r-eps,yc+r+eps];

% generat sample points
N   = 50;
x   = linspace(-pi,pi,N).';
S   = [(r-0.5*eps)*cos(x)+xc,...
       (r-0.5*eps)*sin(x)+yc]; % sample set
   
sdf.Node    = S;
sdf.Element = [(1:N-1).',(2:N).'];

end

function d = dCircle(P,xc,yc,r)
    d = sqrt((P(:,1)-xc).^2+(P(:,2)-yc).^2)-r;
    d=[d,d];
end