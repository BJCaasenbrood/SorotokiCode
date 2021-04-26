function zaxis(Label,varargin)
if isempty(varargin)
    zlabel(Label,'interpreter','latex','fontsize',14);
else
    zlabel([Label, ' ','(',varargin{1},')'],'interpreter','latex','fontsize',14);
end
end