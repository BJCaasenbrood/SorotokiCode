function h = fplot(X,varargin)

if size(X,2) == 2
    h = plot(X(:,1),X(:,2),varargin{:});
elseif size(X,2) == 3 && size(X,3) > 1
    h = plot3(squeeze(X(:,1,:)),squeeze(X(:,2,:)),squeeze(X(:,3,:)),varargin{:});
else
    h = plot3(X(:,1),X(:,2),X(:,3),varargin{:});
end
sorocolor;
end

