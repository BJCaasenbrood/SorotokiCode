function x = Bunny(Request,Arg)

  BdBox = BoundingBox;
  file = 'Bunny';
  
  if nargin < 1, Request = 'Help'; end
  
  switch(Request)
    case('BdBox');   x = BoundingBox;            % call bounding box
    case('Node');    x = LoadNode(file);         % call distance function
    case('Element'); x = LoadElement(file);      % call distance function
    case('Dist');    x = DistanceFunction(Arg);
    case('BC');      x = BoundaryCondition(Arg{:},BdBox);  % call boundary condition
    case('AP');      x = AnchorPoints(Arg{:},BdBox);
    case('PFix');    x = FixedPoints(BdBox);               % call fixed points
    case('Area');    x = (BdBox(2)-BdBox(1))*...           % call bounding area
                         (BdBox(4)-BdBox(3))*...
                         (BdBox(6)-BdBox(5));
    case('Type');    x = 'stl';
    case('Scale');   x = 1;
    case('Dim');     x = 3;
    case('Material');x = 'Porcelain';
    case('View');    x = SetView(240,20);
    case('Check');   x = CheckModel(file);
    case('Help');    x = Help();
    otherwise;       x = [];
  end
  
end 

%-------------------------------------------------- generate bounding box
function BdBox = BoundingBox()
BdBox = [-37, 40, -33, 41, 0, 83];
end  

%----------------------------------------------- compute distance functions
function Dist = DistanceFunction(structure)
FV = structure{1}; P = structure{2};
[Dist] = Point2Mesh(FV, 'QueryPoints', P,...
    'algorithm','linear_vectorized_subfunctions'); 
end

%----------------------------------------------- compute distance functions
function Element = LoadElement(file)
     ModelFile = [file,'Mat.mat'];
     load(ModelFile);
end

%----------------------------------------------- compute distance functions
function x = SetView(a,b)
view(a,b);
x = [];
end

%----------------------------------------------- compute distance functions
function Node = LoadNode(file)
     ModelFile = [file,'Mat.mat'];
     load(ModelFile);
end

% %---------------------------------------------- compute boundary conditions
% function CheckModel(file)
% 
% Name = ['src\example_models\',file,'Mat.mat'];
% 
% if isfile(Name)
% else
%    [Element,Node] = stlreader([file,'.stl']);
%    save([getPath,'\src\examples\models\',file,'Mat.mat'],'Element','Node');
% end
% 
% end

%---------------------------------------------- compute boundary conditions
function x = BoundaryCondition(Node,Element,BdBox)

  Supp = [];
  Load = [];

  x = {Supp,Load};

end

%---------------------------------------------- compute boundary conditions
function x = AnchorPoints(Node,Element,BdBox)

  x1 = 41.52; y1 = 1.885;  z1 = 79.18;
  x2 = 10.58; y2 = -11.87; z2 = 79.54;

  find()
  
  Anchor = [];
  

  x = {Anchor};

end


%----------------------------------------------------- specify fixed points
function PFix = FixedPoints(BdBox)
   PFix = [];
end

%----------------------------------------------------- specify fixed points
function Help
%CallDisplay('commands:');
CallExecuted('BdBox: bounding box');
CallExecuted('PFix: bounding box');
end