function [x,y,z]=rtubeplot(curve,r,rb,s,n,ct,varargin)
% Usage: [x,y,z]=tubeplot(curve,r,n,ct)
% 
% Tubeplot constructs a tube, or warped cylinder, along
% any 3D curve, much like the build in cylinder function.
% If no output are requested, the tube is plotted.
% Otherwise, you can plot by using surf(x,y,z);
%
% Example of use:
% t=linspace(0,2*pi,50);
% tubeplot([cos(t);sin(t);0.2*(t-pi).^2],0.1);
% daspect([1,1,1]); camlight;
%
% Arguments:
% curve: [3,N] vector of curve data
% r      the radius of the tube
% n      number of points to use on circumference. Defaults to 8
% ct     threshold for collapsing points. Defaults to r/2 
%
% The algorithms fails if you have bends beyond 90 degrees.
% Janus H. Wesenberg, july 2004
% Modified 14 Nov, 2022 - Brandon Caasenbrood
% re     the radius at the end
%

  if nargin<3 || isempty(n), n=8;
     if nargin<2, error('Give at least curve and radius');
    end;
  end;
  if size(curve,1)~=3
    error('Malformed curve: should be [3,N]');
  end;
  if nargin<4 || isempty(ct)
    ct=0.5*r;
  end
  
  if isempty(varargin) 
      ramp = 0;
  else
      if ~isempty(varargin{1})
        ramp = varargin{1};
      else
        ramp = 0;
      end
  end
  
  %Collapse points within 0.5 r of each other
  npoints=1;
  for k=2:(size(curve,2)-1)
    if norm(curve(:,k)-curve(:,npoints))>ct
      npoints=npoints+1;
      curve(:,npoints)=curve(:,k);
    end
  end
  %Always include endpoint
  if norm(curve(:,end)-curve(:,npoints))>0
    npoints=npoints+1;
    curve(:,npoints)=curve(:,end);
  end

  %deltavecs: average for internal points.
  %           first strecth for endpoitns.
  dv=curve(:,[2:end,end])-curve(:,[1,1:end-1]);

  %make nvec not parallel to dv(:,1)
  nvec=zeros(3,1);
  [buf,idx]=min(abs(dv(:,1))); nvec(idx)=1;

  xyz=repmat([0],[3,n+1,npoints+2]);
  
  X = linspace(0,2*pi,n+1);
  [Fx,Fy] = squircle(X,s,r,rb);
  
  %precalculate cos and sing factors:
  cfact=repmat(Fx,[3,1]);
  sfact=repmat(Fy,[3,1]);
  
  ramp = min(max(ramp,1e-6),1-1e-6);
  if numel(ramp) == 1
    R = linspace(1,1-ramp,npoints);
  else
    x = linspace(0,1,numel(ramp));  
    y = linspace(0,1,npoints);  
    R = interp1(x,ramp(:),y);  
  end
  
  %Main loop: propagate the normal (nvec) along the tube
  for k=1:npoints
    convec=cross(nvec,dv(:,k));
    convec=convec./norm(convec);
    nvec=cross(dv(:,k),convec);
    nvec=nvec./norm(nvec);
    %update xyz:
    xyz(:,:,k+1)=repmat(curve(:,k),[1,n+1])+...
        cfact.*repmat(R(k)*nvec,[1,n+1])...
        +sfact.*repmat(R(k)*convec,[1,n+1]);
  end
  
  %finally, cap the ends:
  xyz(:,:,1)   =repmat(curve(:,1),[1,n+1]);
  xyz(:,:,end) =repmat(curve(:,end),[1,n+1]);
  
  %,extract results:
  x=squeeze(xyz(1,:,:));
  y=squeeze(xyz(2,:,:));
  z=squeeze(xyz(3,:,:));
end

%%
function [x,y] = squircle(t,s,ra,rb)
%SQUIRCLE squircle(th,s,ra,rb)

% clamp s values
s = min(max(s,1e-6),1-1e-6);

RA = (ra)/(2*s);
RB = (rb)/(2*s);

Rx1 = RA*(2 + 2*s*sqrt(2) * cos(t) + s^2 * cos(2*t)).^0.5;
Ry1 = RB*(2 + 2*s*sqrt(2) * sin(t) - s^2 * cos(2*t)).^0.5;

Rx2 = -RA*(2 - 2*s*sqrt(2) * cos(t) + s^2 * cos(2*t)).^0.5;
Ry2 = -RB*(2 - 2*s*sqrt(2) * sin(t) - s^2 * cos(2*t)).^0.5;

x = Rx1 + Rx2;
y = Ry1 + Ry2;


end


  
  %... and plot:
%   if nargout<3, surf(x,y,z); end;
  