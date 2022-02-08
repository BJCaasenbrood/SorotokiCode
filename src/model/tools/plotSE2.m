function h = plotSE2(g,varargin)
if isempty(varargin)
   ID = [1 3];  % x-z plane 
elseif strcmp(varargin{1},'xz') || strcmp(varargin{1},'zx')
   ID = [1 3];  % x-z plane 
elseif strcmp(varargin{1},'xy') || strcmp(varargin{1},'yx')
   ID = [1 2];  % x-z plane 
elseif strcmp(varargin{1},'zy') || strcmp(varargin{1},'yz')
   ID = [2 3];  % x-z plane 
end

I = eye(3);
I = I(ID,:);
singleFrame = true;

if numel(g) == 7
    R0 = quat2rot(g(1:4)); 
    p0 = g(5:7);
    R = I*R0;
    p = (p0.'*I')';
elseif numel(g) == 16
    R0 = g(1:3,1:3)*I';
    p0 = g(1:3,4);    
    R = R0(ID,:);
    p = (p0.'*I')';    
elseif numel(g) == 9
    R = g(1:2,1:2);
    p = g(1:2,3);       
elseif isflint(numel(g)/16)
    p0  = reshape(g(1:3,4,:),3,[]);
    Rx0 = reshape(g(1:3,1,:),3,[]);
    Ry0 = reshape(g(1:3,2,:),3,[]);
    Rz0 = reshape(g(1:3,3,:),3,[]);
    p = (p0.'*I');
    Rx = (Rx0.'*I');
    Ry = (Ry0.'*I');
    Rz = (Rz0.'*I');
    
    singleFrame = false;
end

if size(varargin,2) < 2
    a = 1;
else
    a = varargin{2};
end
h = cell(2,1);
if singleFrame
    for ii = 1:2
        h{ii} = plotarrow(p,a*R(:,ii),ii);
    end
else
   for ii = 1:2
        s = ID(ii);
        if s == 1,     h{ii} = plotarrowmany(p,a*Rx,ii);
        elseif s == 2, h{ii} = plotarrowmany(p,a*Ry,ii);  
        else,          h{ii} = plotarrowmany(p,a*Rz,ii);    
        end
   end
end

end

function h = plotarrow(p,n,id)
hold on; 
h = quiver(p(1),p(2),...
    n(1), n(2),'Color',col(id));
end

function h = plotarrowmany(p,n,id)
hold on; 
h = quiver(p(:,1),p(:,2),...
    n(:,1), n(:,2),'Color',col(id));
end


