function h = fplot(X,varargin)

if size(X,2) == 2
    h = plot(X(:,1),X(:,2),varargin{:});
else
    h = plot3(X(:,1),X(:,2),X(:,3),varargin{:});
end
sorocolor;
end

