function SetQuiverLength(q,mags,varargin)
%--------------------------------------------------
% function SetQuiverLength(q,mags)
%
%   Function that sets the length of the quiver
%   You can give a list of length (in x,y axis units), the function will rescale the vectors to the right length
%
% Input: q    = handle to quiver plot, can be quiver or quiver3
%        mags = Desired length of each vector in units of the x,y axis
%               If mags is a scalar, all the vector will have that length
% Optional inputs: given as ('variable',value)
%        'HeadLength' = single value for the length of the head in units of the xyz axis.
%                       It is applied to all the vectors
%        'HeadAngle' = angle between the two lines forming the head
%                      default=28.0724^\circ
%        'RotHead'   = Angle [deg] by which the head will be rotated around the vector axis.
%                      This allows to set the head in different planes
%
% NOTE: For some unknown reason, MATLAB does not always simply change the length of the vectors and requires a
%       pause statement towards the end.
%       Would your vectors not be the right size, increase the duration of the pause
%       (or even better, suggest a solution that does not require a pause)
%
% Example:
%     [x,y] = meshgrid(-2:.2:2,-1:.15:1);
%     z = x .* exp(-x.^2 - y.^2);
%     [u,v,w] = surfnorm(x,y,z);
%     q = quiver3(x,y,z,u,v,w); hold on; surf(x,y,z); hold off;
%     drawnow;                  % This is needed as, if the plot is not plotted, the VertexData for the quiver are not existent.
%     view(180,0);
%     mag = 0.1*ones(size(u));  % Length of the vectors : 0.1
%     SetQuiverLength(q,mag,'HeadLength',0.05,'HeadAngle',90);
%
%--------------------------------------------------
%   Authorship:
%     Function made by Alexandre De Spiegeleer. Feel free to use it.
%--------------------------------------------------

%// Set default values of varargin
HeadLength = [];
HeadAngle = [];
RotHead = [];

%// Read the optional inputs
if find(strcmp('HeadLength',varargin))
  HeadLength = varargin{ find(strcmp('HeadLength',varargin))+1 };
end
if find(strcmp('HeadAngle',varargin))
  HeadAngle = varargin{ find(strcmp('HeadAngle',varargin))+1 };
end
if find(strcmp('RotHead',varargin))
  RotHead = varargin{ find(strcmp('RotHead',varargin))+1 };
end

%// Start by removing the autoscale option
set(q,'AutoScale','Off');

%// Get the Head and Tail options
H = get(q.Head);
T = get(q.Tail);

%// Handle mags input
if isempty(mags)
  warning('No input was given for the length of the vectors. They are now all of length 1.');
  mags = ones(size(T.VertexData,2)/2,1);
elseif numel(mags) == 1
  mags = repmat(mags,size(T.VertexData,2)/2,1);
elseif size(mags) == size(q.XData)
  mags = reshape(mags,[],1);
elseif isrow(mags)
  mags = mags';
elseif numel(mags)~=size(T.VertexData,2)/2
  warning('Different number of magnitudes and quiver vectors!');
  return;
else
  warning('Bug when calling SetQuiverLength function!');
  return;
end

if find(mags<0)
  warning(['At least one of the length is requested to be negative. ' ...
          'This will swap the direction of the vector(s)!']);
end

%// The data for which there is no vector should be removed
%// as there is no Head or Tail data for those.
if numel(mags) == numel(q.UData)
  mags(isnan(q.UData)) = [];
end

%// Reshape the head and the tail (It is easier to think about it then but requires some repmat and permute)
%// The head has 3 vertices and the tail has 2.
Tail_ori = reshape(T.VertexData,size(T.VertexData,1),2,[]);
Head_ori = reshape(H.VertexData,size(H.VertexData,1),3,[]);

%// Measure original length of Head and Tail
Tail_length = squeeze(sqrt(sum(diff(Tail_ori,1,2).^2,1)));
Head_length = squeeze(sqrt(sum(diff(Head_ori(:,[1,2],:),1,2).^2,1)));

%// Head-Tail original ratio
Length_ratio = Head_length./Tail_length;

%// New length of head
if isempty(HeadLength)
  Head_mag = mags.*Length_ratio;
  Head_mag = permute(Head_mag(:,:,ones(size(H.VertexData,1),1)),[3 2 1]);
elseif ~isempty(HeadLength)
  Head_mag = repmat(HeadLength,[size(H.VertexData,1) 1 size(Head_ori,3)]);
end

%// Direction tail
Tail_dir = diff(Tail_ori,1,2)./permute(Tail_length(:,:,ones(size(T.VertexData,1),1)),[3 2 1]);

%// Directions of head
Head_dir_a = diff(Head_ori(:,[2,1],:),1,2)./permute(Head_length(:,:,ones(size(H.VertexData,1),1)),[3 2 1]);
Head_dir_b = diff(Head_ori(:,[2,3],:),1,2)./permute(Head_length(:,:,ones(size(H.VertexData,1),1)),[3 2 1]);

%// Set new HeadAngle by use of Rodrigues' rotation formula
if ~isempty(HeadAngle)
  % rotation axis:
  k = cross(Head_dir_b,Head_dir_a,1)./sqrt(sum(cross(Head_dir_b,Head_dir_a,1).^2,1));
  % Original angle:
  alpha = acosd( dot(Head_dir_a,Head_dir_b,1) );
  % Angle of rotation:
  theta = (HeadAngle-alpha)/2;
  % Rotation
  ct = repmat(cosd(theta),size(H.VertexData,1),1,1);
  st = repmat(sind(theta),size(H.VertexData,1),1,1);
  Head_dir_a_r = Head_dir_a .* ct + ...
                cross(k,Head_dir_a,1) .* st + ...
                k .* repmat(dot(k,Head_dir_a,1),size(H.VertexData,1),1,1) .* (1-ct);
  ct = repmat(cosd(-theta),size(H.VertexData,1),1,1);
  st = repmat(sind(-theta),size(H.VertexData,1),1,1);
  Head_dir_b_r = Head_dir_b .* ct + ...
                cross(k,Head_dir_b,1) .* st + ...
                k .* repmat(dot(k,Head_dir_b,1),size(H.VertexData,1),1,1) .* (1-ct);
  Head_dir_a = Head_dir_a_r;
  Head_dir_b = Head_dir_b_r;
end

%// Set new RotHead by use of Rodrigues' rotation formula
if ~isempty(RotHead)
  % rotation axis:
  k = Tail_dir;
  % Angle of rotation:
  theta = RotHead;
  % Rotation
  ct = repmat(cosd(theta),size(H.VertexData,1),1,1);
  st = repmat(sind(theta),size(H.VertexData,1),1,1);
  Head_dir_a_r = Head_dir_a .* ct + ...
                cross(k,Head_dir_a,1) .* st + ...
                k .* repmat(dot(k,Head_dir_a,1),size(H.VertexData,1),1,1) .* (1-ct);
  Head_dir_b_r = Head_dir_b .* ct + ...
                cross(k,Head_dir_b,1) .* st + ...
                k .* repmat(dot(k,Head_dir_b,1),size(H.VertexData,1),1,1) .* (1-ct);
  Head_dir_a = Head_dir_a_r;
  Head_dir_b = Head_dir_b_r;
end

%// New Tail End = start of tail + magnitude*normalized direction
TailVertex = Tail_ori;
TailVertex(:,2,:) = TailVertex(:,1,:) +  ...
            permute(mags(:,:,ones(size(T.VertexData,1),1)),[3 2 1]) .*  Tail_dir;

%// Middle point of Head is same as end of Tail
HeadVertex = Head_ori;
HeadVertex(:,2,:) = TailVertex(:,2,:);

%// Head vector a
HeadVertex(:,1,:) = HeadVertex(:,2,:) + Head_mag .* Head_dir_a;
%// Head vector b
HeadVertex(:,3,:) = HeadVertex(:,2,:) + Head_mag .* Head_dir_b;


%// If the vector is originally of size 0, the vectors are now NaN.
%// Change those back to have just a simple "single point"
Idx_NaN = find(isnan(squeeze(TailVertex(1,end,:))));
TailVertex(:,end,Idx_NaN) = TailVertex(:,1,Idx_NaN);
HeadVertex(:,:,Idx_NaN) = repmat(TailVertex(:,1,Idx_NaN),1,3,1);

%// Reshape the data
HeadVertex = reshape(HeadVertex,size(H.VertexData,1),[]);
TailVertex = reshape(TailVertex,size(T.VertexData,1),[]);

%// Don't ask me about this line ...
%// Matlab sucks and does not do the change graphically if I don't have a pause or a keyboard ...
pause(0.0102)

%// Set the data to the quiver properties
set(q.Head,'VertexData',HeadVertex);
set(q.Tail,'VertexData',TailVertex);
