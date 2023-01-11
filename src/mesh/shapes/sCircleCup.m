function sdf = sCircleCup(xc,yc,r,x)

if nargin < 2
    r = xc; 
    xc = 0;
    yc = 0;
    x = r;
end

eps = 1e-4*r;

if x > 0
    sdf = Sdf(@(P) dDiff(dRectangle(P,xc-r,xc+r,yc-x,yc),...
        dCircle(P,xc,yc,r)));

    sdf.BdBox = [xc-r-eps,xc+r+eps,yc-x-eps,yc+r+eps];
else
    sdf = Sdf(@(P) dUnion(dRectangle(P,xc-r,xc+r,yc,yc+abs(x)),...
        dCircle(P,xc,yc,r)));

    sdf.BdBox = [xc-r-eps,xc+r+eps,yc-r-eps,yc+abs(x)+eps];
end

end

function d = dCircle(P,xc,yc,r)
    d = sqrt((P(:,1)-xc).^2+(P(:,2)-yc).^2)-r;
    %d = [d,d];
end

function d = dRectangle(P,x1,x2,y1,y2)
d = [x1-P(:,1), P(:,1)-x2, y1-P(:,2), P(:,2)-y2];
d = max(d,[],2);
end