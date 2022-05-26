function h = fquiver(X,N,varargin)

if size(X,2) == 2
    h = quiver(X(:,1),X(:,2),...
        N(:,1),N(:,2),varargin{:});
else
    h = quiver3(X(:,1),X(:,2),X(:,3),...
        N(:,1),N(:,2),N(:,3),varargin{:});
end
end
