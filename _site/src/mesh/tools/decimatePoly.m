function [C_out,i_rem,CI]=decimatePoly(C,opt,vis)
% Reduce complexity of a 2D simple (i.e. non-self intersecting), closed 
% piecewise linear contour by specifying boundary offset tolerance.
%
% ******************************* IMPORTANT *******************************
% This function does not guarantee that simplified contour will be free of
% self-intersections.
% *************************************************************************
%
%
% INPUT:
%   - C     : N-by-2 array containing an ordered list of vertex 
%             co-ordinates, such that the first, C(1,:), and last, 
%             C(end,:), vertices are the same. 
%   - opt   : opt can be specified in one of two ways: 
%             ----------------------APPROACH #1 {default} -----------------
%             - opt : opt=[B_tol 1], where B_tol is the maximum acceptable 
%                     offset from the original boundary, B_tol must be 
%                     expressed in the same length units as the co-ords in 
%                     C. Default setting is B_tol=Emin/2, where Emin is the
%                     length of the shortest edge.  
%             ----------------------APPROACH #2----------------------------
%              - opt : opt=[P_tol 2], where P_tol is the fraction of the 
%                      total number of contour vertices to be retained. 
%                      Accordingly, P_tol must be a real number on the
%                      interval (0,1).
%   - vis   : binary option used to specify whether to print out the effect
%             of decimation on contour area and perimeter. vis=true is the 
%             default setting. Set vis=false to suppress the print out.
%
% OUTPUT:
%   - C_out : M-by-2 array containing co-ordinates of simplified contour.
%   - i_rem : N-by-1 logical array used to mark vertices of the original 
%             contour that were removed by decimation.
%   - CI    : structure containing summary "parameters"  (i.e., # of 
%             vertices, perimeter, and area) of the input and output
%             contours.
%
% ALGORITHM:
% (1) For every vertex compute boundary offset error.  
% (2) Rank all vertices in increasing order of error.
% (3) Remove vertex with the lowest error.
% (4) Recompute and accumulate errors for two vertices adjacent to the 
%     deleted vertex and go back to (2). 
% (5) Repeat (2) to (4) until no further vertices can be removed or the 
%     number of vertices has reached a user-defined value.
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%
% Check input args
if nargin<2, opt={}; end
if nargin<3 || isempty(vis), vis=true; end
opt=CheckInputArgs(C,opt,vis);
N=size(C,1);
i_rem=false(N,1); 
if N<=4, 
    C_out=C; 
    return
end
% Tolerance parameter, perimeter and area of the input contour
[Po,Emin]=PolyPerim(C);
B_tol=Emin/2;
Ao=PolyArea(C);
No=N-1;
if ~isempty(opt), B_tol=opt(1); end
if isempty(opt), opt=[B_tol 1]; end
Nmin=3;
if opt(2)==2
    Nmin=round((N-1)*opt(1));
    if (N-1)==Nmin, return; end
    if Nmin<3, Nmin=3; end
end
% Remove (repeating) end-point
C(end,:)=[];
N=N-1;
% Compute distance offset errors ------------------------------------------
D31=circshift(C,[-1 0])-circshift(C,[1 0]);
D21=C-circshift(C,[1 0]);
dE_new2=sum(D31.^2,2); % length^2 of potential new edges
% Find closest points to current vertices on the new edges
t=sum(D21.*D31,2)./dE_new2; 
t(t<0)=0;
t(t>1)=1;
V=circshift(C,[1 0])+bsxfun(@times,t,D31);
% Evaluate distance^2
Err_D2=sum((V-C).^2,2);
% Initialize distance error accumulation array 
DEAA=zeros(N,1);
% Begin decimation --------------------------------------------------------
idx_ret=1:N; % keep track of retained vertices
while true
    
    % Find vertices whose removal will satisfy the decimation criterion
    idx_i=Err_D2<B_tol;
    if sum(idx_i)==0 && N>Nmin && opt(2)==2
        B_tol=B_tol*sqrt(1.5);
        continue
    end
    idx_i=find(idx_i);
    if isempty(idx_i) || N==Nmin, break; end
    N=N-1;
    
    % Vertex with the smallest net error
    [~,i_min]=min(Err_D2(idx_i));
    idx_i=idx_i(i_min);
    
    % Update distance error accumulation array 
    DEAA(idx_i)=DEAA(idx_i)+sqrt(Err_D2(idx_i));
    
    i1=idx_i-1; if i1<1, i1=N; end
    i3=idx_i+1; if i3>N, i3=1; end
    
    DEAA(i1)=DEAA(idx_i);
    DEAA(i3)=DEAA(idx_i);
    
    % Recompute errors of vertices adjacent to vertex marked for deletion
    i1_1=i1-1; if i1_1<1, i1_1=N; end
    i1_3=i3;
    
    i3_1=i1;
    i3_3=i3+1; if i3_3>N, i3_3=1; end
    
    err_D1=RecomputeErrors(C([i1_1,i1,i1_3],:));
    err_D3=RecomputeErrors(C([i3_1,i3,i3_3],:));
    
    % Upadate errors
    Err_D2(i1)=(sqrt(err_D1)+ DEAA(i1)).^2;
    Err_D2(i3)=(sqrt(err_D3)+ DEAA(i3)).^2;
        
    % Remove vertex
    C(idx_i,:)=[];
    idx_ret(idx_i)=[];
    DEAA(idx_i)=[];
    
    Err_D2(idx_i)=[];
    
end
C=[C;C(1,:)]; C_out=C;
i_rem(idx_ret)=true;
i_rem=~i_rem;
i_rem(end)=i_rem(1);
% Perimeter and area of simplified contour
P=PolyPerim(C);
A=PolyArea(C);
% Contour info
CI.in.num_verts=No;
CI.in.perimeter=Po;
CI.in.area=Ao;
CI.out.num_verts=N;
CI.out.perimeter=P;
CI.out.area=A;
% Effect of decimation on total area and perimeter of the contour
if ~vis, return; end
fprintf('%+20s%+15s%+11s\n','# of verts','perimeter','area')
fprintf('-----------------------------------------------------\n')
fprintf('%-10s%-16u%-16.3f%-16.3f\n','in',No,Po,Ao)
fprintf('%-10s%-16u%-16.3f%-16.3f\n','out',N,P,A)
fprintf('-----------------------------------------------------\n')
fprintf('%-10s%-16.3f%-16.3f%-16.3f\n\n','% change',(N-No)/No*100,(P-Po)/Po*100,(A-Ao)/Ao*100)
%==========================================================================
function err_D2=RecomputeErrors(V)
% Recompute distance offset error for a small subset of vertices.  
%
%   - V     : 3-by-2 array of triangle vertices, where V(2,:) is the vertex
%             marked for removal.
% Distance offset error
D31=V(3,:)-V(1,:);
D21=V(2,:)-V(1,:);
dE_new2=sum(D31.^2,2); % length^2 of potential new edge
% Find closest point to the current vertex on the new edge
t=sum(D21.*D31,2)/dE_new2; 
t(t<0)=0;
t(t>1)=1;
p=V(1,:)+bsxfun(@times,t,D31);
% Evaluate the distance^2
err_D2=sum((p-V(2,:)).^2);
%==========================================================================
function A=PolyArea(C)
% Polygon area
dx=C(1:(end-1),1)-C(2:end,1);
dy=C(1:(end-1),2)+C(2:end,2);
A=sum(dx.*dy)/2;
A=abs(A);
%==========================================================================
function [P,Emin]=PolyPerim(C)
% Polygon perimeter
dE=C(2:end,:)-C(1:(end-1),:);
dE=sqrt(sum(dE.^2,2));
P=sum(dE);
Emin=min(dE);
%==========================================================================
function opt=CheckInputArgs(C,opt,vis)
% Check validity of the input arguments
siz=size(C);
if numel(siz)~=2 || siz(2)>siz(1) || ~isnumeric(C) || ~ismatrix(C)
    error('1st input argument (C) must be a N-by-2 array')
end
if norm(C(1,:)-C(end,:))>1E-12
    error('First and last points in C must be the same')
end
if isempty(opt), return; end
if ~isnumeric(opt) || numel(opt)~=2
    error('Invalid entry for 2nd input argument (opt)')
end
if ~(opt(2)==1 || opt(2)==2)
    error('Invalid entry for 2nd input argument (opt). opt(2) must be set to 1 or 2.')
end
if opt(2)==1 && opt(1)<=eps
    error('Invalid entry for 2nd input argument (opt). When opt(2)==1, opt(1) must be greater than zero.')
end
if opt(2)==2 && (opt(1)<=eps || opt(1)>=1)
    error('Invalid entry for 2nd input argument (opt). When opt(2)==2, opt(1) must be on the interval (0,1).')
end
if numel(vis)~=1 || ~islogical(vis)
    error('Invalid entry for 3rd input argument (vis). vis must be binary.')
end
