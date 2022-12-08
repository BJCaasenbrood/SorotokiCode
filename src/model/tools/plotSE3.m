function h = plotSE3(H,varargin)

arrowsize = 1;
linewidth = 2;

if size(H,2) == 7 || size(H,1) == 7
    p = H(5:7);
    R = quat2rot(H(1:4));
else
    p = H(1:3,4);
    R = H(1:3,1:3);
end

if ~isempty(varargin)
    for ii = 2:2:numel(varargin)
        if strcmp(varargin{ii-1},'linewidth')
            linewidth = varargin{ii};
        elseif strcmp(varargin{ii-1},'arrowsize')
            arrowsize = varargin{ii};
        end
    end
end

% if nargin < 2
%    a = 1; 
% else
%    a = varargin{1};
% end

for ii = 1:3
    hold on; 
    n = arrowsize*R(:,ii);
    h = quiver3(p(1),p(2),p(3),...
    n(1), n(2), n(3),'Color',col(ii),'LineWidth',linewidth);
    
end

end

