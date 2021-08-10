function [d_min, varargout] = p_poly_dist(Pp, Pv, varargin)
%p_poly_dist Find minimum distances from points to a polyline or to a 
% closed polygon.
%
% Description:
% Compute the distances from each one of a set of np points p(1), p(2),...
% p(np) on a 2D plane to a polyline or a closed polygon. Polyline is 
% defined as a set of nv-1 segments connecting nv ordered vertices v(1), 
% v(2), ..., v(nv). The polyline can optionally be treated as a closed
% polygon.
% Distance from point j to a segment k is defined as a distance from this
% point to a straight line passing through vertices v(k) and v(k+1), when
% the projection of point j on this line falls INSIDE segment k; and to
% the closest of v(k) or v(k+1) vertices, when the projection falls OUTSIDE
% segment k. Distance from point j to a polyline is defined as a minimum of
% this point's distances to all segments.
% In a case when the projected point fall OUTSIDE of all polyline segments,
% the distance to a closest vertex of the polyline is returned
%
% Input arguments:
% Required:
% [d_min, varargout] = p_poly_dist(xp, yp, xv, yv)
% xp - vector of points X coordinates (1 X np)
% yp - vector of points Y coordinates (1 X np)
% xv - vector of polygon vertices' X coordinates (1 X nv)
% yv - vector of polygon vertices' Y coordinates (1 X nv)
%
% Optional:
% [d_min, varargout] = p_poly_dist(xp, yp, xv, yv, find_in_out)
%
% find_in_out - logical flag. When true, the polyline is treated as a
% closed polygon, and each input point is checked wether it is inside or 
% outside of this polygon. In such case, an open polyline is automatically 
% closed by adding a last point identical to a first one. 
% Note: when this function is called with ver. 1.0 signature, i.e:
% d = p_poly_dist(xp, yp, xv, yv)
% the flag find_in_out is assumed to be true, to keep the functionality
% compatible with a previous version.
%
% Output arguments:
% Required:
% d_min = p_poly_dist(xp, yp, xv, yv, varargin)
%
% d_min - vector of distances (1 X np) from points to polyline. This is 
% defined as either a distance from a point to it's projection on a line 
% that passes through a pair of consecutive polyline vertices, when this 
% projection falls inside a segment; or as a distance from a point to a 
% closest vertex, when the projection falls outside of a segment. When 
% find_in_out input flag is true, the polyline is assumed to be a closed 
% polygon, and distances of points INSIDE this polygon are defined as 
% negative.
% 
% Optional:
% [d_min, x_d_min, y_d_min, is_vertex, xc, yc, idx_c, Cer, Ppr] = 
%       p_poly_dist(xp, yp, xv, yv, varargin);
%
% x_d_min - vector (1 X np) of X coordinates of the closest points of
% polyline. 
%
% y_d_min - vector (1 X np) of Y coordinates of the closest points of
% polyline. 
%
% is_vertex - logical vector (1 X np). If is_vertex(j)==true, the closest
% polyline point to a point (xp(j), yp(j)) is a vertex. Otherwise, this
% point is a projection on polyline's segment.
%
% idx_c - vector (1 X np) of indices of segments that contain the closest
% point. For instance,  idx_c(2) == 4 means that the polyline point closest
% to point 2 belongs to segment 4
%
% xc - an array (np X nv-1) containing X coordinates of all projected
% points. xc(j,k) is an X coordinate of a projection of point j on a
% segment k
%
% yc - an array (np X nv-1) containing Y coordinates of all projected
% points. yc(j,k) is Y coordinate of a projection of point j on a
% segment k
%
% is_in_seg - logical array (np X nv-1). If is_in_seg(j,k) == true, the
% projected point j with coordinates (xc(j,k), yc(j,k)) lies INSIDE the
% segment k
%
% Cer - a 3D array (2 X 2 X nv-1). Each 2 X 2 slice represents a rotation
% matrix from an input Cartesian coordinate system to a system defined
% by polyline segments.
%
% Ppr - 3D array of size 2 X np X (nv-1). Ppr(1,j,k) is an X coordinate
% of point j in coordinate systems defined by segment k. Ppr(2,j,k) is its
% Y coordinate.
%
% Routines: p_poly_dist.m
% Revision history:
% Oct 2, 2015 - version 2.0 (Michael Yoshpe). The original ver. 1.0 
% function was completely redesigned. The main changes introduced in 
% ver. 2.0:
% 1. Distances to polyline (instead of a closed polygon in ver. 1.0) are 
% returned. The polyline can optionally be treated as a closed polygon.
% 2. Distances from multiple points to the same polyline can be found
% 3. The algorithm for finding the closest points is now based on 
% coordinate system transformation. The new algorithm avoids numerical 
% problems that ver. 1.0 algorithm could experience in "ill-conditioned" 
% cases.
% 4. Many optional output variables were added. In particular, the closest
% points on polyline can be returned.
% 5. Added input validity checks
% 7/9/2006  - case when all projections are outside of polygon ribs
% 23/5/2004 - created by Michael Yoshpe 
% Remarks:
%**************************************************************************

xv = Pv(:,1); yv = Pv(:,2); 
xp = Pp(:,1); yp = Pp(:,2); 

find_in_out = false;

if(nargin >= 5)
   find_in_out = varargin{1};  
elseif((nargin==4) && (nargout==1)) % mimic ver. 1.0 functionality
   find_in_out = true;
end

% number of points and number of vertices in polyline
nv = length(xv);
np = length(xp);

if(nv < 2)
   error('Polyline must have at least 2 vertices');
end

if((find_in_out == true) && (nv < 3))
   error('Polygon must have at least 3 vertices');
end

% if finding wether the points are inside or outsite the polygon is
% required, make sure the verices represent a closed polygon
if(find_in_out)
   % If (xv,yv) is not closed, close it.
   nv = length(xv);
   if ((xv(1) ~= xv(nv)) || (yv(1) ~= yv(nv)))
      xv = [xv(:)' xv(1)];
      yv = [yv(:)' yv(1)];
      nv = nv + 1;
   end
end

% Cartesian coordinates of the polyline vertices
Pv = [xv(:) yv(:)];

% Cartesian coordinates of the given points
Pp = [xp(:) yp(:)];

% matrix of distances between all points and all vertices
% dpv(j,k) - distance from point j to vertex k
dpv(:,:) = hypot((repmat(Pv(:,1)', [np 1])-repmat(Pp(:,1), [1 nv])),...
                 (repmat(Pv(:,2)', [np 1])-repmat(Pp(:,2), [1 nv])));

% Find the vector of minimum distances to vertices. 
[dpv_min, I_dpv_min] = min(abs(dpv),[],2);

% coordinates of consecutive vertices
P1 = Pv(1:(end-1),:);
P2 = Pv(2:end,:);
dv = P2 - P1;

% vector of distances between each pair of consecutive vertices
vds = hypot(dv(:,1), dv(:,2));

% check for identical points
idx = find(vds < 10*eps);
if(~isempty(idx))
   error(['Points ' num2str(idx) ' of the polyline are identical']);
end

% check for a case when closed polygon's vertices lie on a stright line, 
% i.e. the distance between the last and first vertices is equal to the sum
% of all segments except the last
if(find_in_out)
   s = cumsum(vds);
   if((s(end-1) - vds(end)) < 10*eps)
      error('Polygon vertices should not lie on a straight line');
   end
end

% Each pair of consecutive vertices P1(j), P2(j) defines a rotated 
% coordinate system with origin at P1(j), and x axis along the vector 
% (P2(j)-P1(j)). 
% Build the rotation matrix Cer from original to rotated system
ctheta = dv(:,1)./vds;
stheta = dv(:,2)./vds;
Cer = zeros(2,2,nv-1);
Cer(1,1,:) = ctheta;
Cer(1,2,:) = stheta;
Cer(2,1,:) = -stheta;
Cer(2,2,:) = ctheta;

% Build the origin translation vector P1r in rotated frame by rotating the
% P1 vector
P1r = [(ctheta.*P1(:,1) + stheta.*P1(:,2)),...
       -stheta.*P1(:,1) + ctheta.*P1(:,2)];

Cer21 = zeros(2, nv-1);
Cer22 = zeros(2, nv-1);

Cer21(:,:) = Cer(1,:,:);
Cer22(:,:) = Cer(2,:,:);

% Ppr is a 3D array of size 2 * np * (nv-1). Ppr(1,j,k) is an X coordinate
% of point j in coordinate systems defined by segment k. Ppr(2,j,k) is its
% Y coordinate.

% Rotation. Doing it manually, since Matlab cannot multiply 2D slices of ND
% arrays.
Ppr(1,:,:) = Pp*Cer21;
Ppr(2,:,:) = Pp*Cer22;

% translation
Ppr(1,:,:) = Ppr(1,:,:) - permute(repmat(P1r(:,1), [1 1 np]), [2 3 1]);
Ppr(2,:,:) = Ppr(2,:,:) - permute(repmat(P1r(:,2), [1 1 np]), [2 3 1]);

% Pcr is a 3D array of size 2 * np * (nv-1) that holds the projections of
% points on X axis of rotated coordinate systems. Pcr(1,j,k) is an X
% coordinate of point j in coordinate systems defined by segment k.
% Pcr(2,j,k) is its Y coordinate, which is identically zero for projected
% point.
Pcr = zeros(size(Ppr));
Pcr(1, :, :) = Ppr(1,:,:);
Pcr(2, :, :) = 0;

% Cre is a rotation matrix from rotated to original system
Cre = permute(Cer, [2 1 3]);

% Pce is a 3D array of size 2 * np * (nv-1) that holds the projections of
% points on a segment in original coordinate systems. Pce(1,j,k) is an X
% coordinate of the projection of point j on segment k.
% Pce(2,j,k) is its Y coordinate
Pce = zeros(2,np,(nv-1));
Pce(1,:,:) = Pcr(1,:,:).*repmat(Cre(1,1,:), [1 np 1]) + ...
             Pcr(2,:,:).*repmat(Cre(1,2,:), [1 np 1]);
Pce(2,:,:) = Pcr(1,:,:).*repmat(Cre(2,1,:), [1 np 1]) + ...
             Pcr(2,:,:).*repmat(Cre(2,2,:), [1 np 1]);
          
% Adding the P1 vector
Pce(1,:,:) = Pce(1,:,:) + permute(repmat(P1(:,1), [1 1 np]), [2 3 1]);
Pce(2,:,:) = Pce(2,:,:) + permute(repmat(P1(:,2), [1 1 np]), [2 3 1]);

% x and y coordinates of the projected (cross-over) points in original
% coordinate frame
xc = zeros(np, (nv-1));
yc = zeros(np, (nv-1));
xc(:,:) = Pce(1,:,:);
yc(:,:) = Pce(2,:,:);

r = zeros(np,(nv-1));
cr = zeros(np,(nv-1));
r(:,:) = Ppr(1,:,:);
cr(:,:) = Ppr(2,:,:);

% Find the projections that fall inside the segments
is_in_seg = (r>0) & (r<repmat(vds(:)', [np  1 1]));

% find the minimum distances from points to their projections that fall
% inside the segments (note, that for some points these might not exist,
% which means that closest points are vertices)
B = NaN(np,nv-1);
B(is_in_seg) = cr(is_in_seg);

[cr_min, I_cr_min] = min(abs(B),[],2);

% Build the logical index which is true for members that are vertices,
% false otherwise. Case of NaN in cr_min means that projected point falls
% outside of a segment, so in such case the closest point is vertex.

% point's projection falls outside of ALL segments
cond1 = isnan(cr_min); 

% point's projection falls inside some segments, find segments for which
% the minimum distance to vertex is smaller than to projection
cond2 = ((I_cr_min ~= I_dpv_min) & (cr_min - dpv_min)> 0);

is_vertex = (cond1 | cond2);

% build the minimum distances vector
d_min = cr_min;
d_min(is_vertex) = dpv_min(is_vertex);

% mimic the functionality of ver. 1.0 - make all distances negative for
% points INSIDE the polygon
if(find_in_out)
   in = inpolygon(xp, yp, xv, yv);
   d_min(in) = -d_min(in);
end

% initialize the vectors of minimum distances to the closest vertices
nz = max(np, nv);

vtmp = zeros(nz, 1);
vtmp(1:nv) = xv(:);
x_d_min = vtmp(I_dpv_min);

vtmp = zeros(nz, 1);
vtmp(1:nv) = yv(:);
y_d_min = vtmp(I_dpv_min);

% replace the minimum distances with those to projected points that fall
% inside the segments
idx_pr = sub2ind(size(xc), find(~is_vertex), I_cr_min(~is_vertex));
x_d_min(~is_vertex) = xc(idx_pr);
y_d_min(~is_vertex) = yc(idx_pr);

% find the indices of segments that contain the closest points
% note that I_dpv_min contains indices of closest POINTS. To find the 
% indices of closest SEGMENTS, we have to substract 1
idx_c = I_dpv_min-1; 
[ii,jj] = ind2sub(size(xc), idx_pr);
idx_c(ii) = jj;

% assign optional outputs
switch nargout
   case 0
   case 1
   case 2
      varargout{1} = x_d_min;
   case 3
      varargout{1} = x_d_min;
      varargout{2} = y_d_min;
   case 4
      varargout{1} = x_d_min;
      varargout{2} = y_d_min;
      varargout{3} = is_vertex;
   case 5
      varargout{1} = x_d_min;
      varargout{2} = y_d_min;
      varargout{3} = is_vertex;
      varargout{4} = idx_c;
   case 6
      varargout{1} = x_d_min;
      varargout{2} = y_d_min;
      varargout{3} = is_vertex;
      varargout{4} = idx_c;
      varargout{5} = xc;
   case 7
      varargout{1} = x_d_min;
      varargout{2} = y_d_min;
      varargout{3} = is_vertex;
      varargout{4} = idx_c;
      varargout{5} = xc;
      varargout{6} = yc;
   case 8
      varargout{1} = x_d_min;
      varargout{2} = y_d_min;
      varargout{3} = is_vertex;
      varargout{4} = idx_c;
      varargout{5} = xc;
      varargout{6} = yc;
      varargout{7} = is_in_seg;
   case 9
      varargout{1} = x_d_min;
      varargout{2} = y_d_min;
      varargout{3} = is_vertex;
      varargout{4} = idx_c;
      varargout{5} = xc;
      varargout{6} = yc;
      varargout{7} = is_in_seg;      
      varargout{8} = Cer;
   case 10
      varargout{1} = x_d_min;
      varargout{2} = y_d_min;
      varargout{3} = is_vertex;
      varargout{4} = idx_c;
      varargout{5} = xc;
      varargout{6} = yc;
      varargout{7} = is_in_seg;            
      varargout{8} = Cer;
      varargout{9} = Ppr;      
   otherwise
      error('Invalid number of output arguments, must be between 1 and 9');
end

end