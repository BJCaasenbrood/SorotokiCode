function P = DijkstraFilter(Fem)

tri = ConcaveDelaunay(Fem,0.95);

A = zeros(Fem.NElem);
I = tri(:); 
J = tri(:,[2 3 1]); J = J(:);
IJ = I + Fem.NElem*(J-1); 
A(IJ) = 1;

[Dist,~] = Dijkstra(A,Fem.Center);

P = spdiags(1./max(Dist,[],1).',0,Fem.NElem,Fem.NElem)*Dist;

end

function tri = ConcaveDelaunay(fem,per)

Node = fem.Node;
X = fem.Center(:,1);
Y = fem.Center(:,2);
xy = fem.Center;
tri = delaunay(X,Y);

%maxEdge = 0.1;
edges = [tri(:,[1 2]);tri(:,[1 3]);tri(:,[2 3])]; 
CellEdge = num2cell(edges,2);

NodePos = @(E) sqrt((Node(E(1),1)-Node(E(2),1))^2 + ...
    (Node(E(1),1)-Node(E(2),1))^2 );  

CellDc = cellfun(NodePos,CellEdge,'UniformOutput',false); 
dC = sort(vertcat(CellDc{:})); dCm = mean(dC(1:100));

maxEdge = 10*dCm;

 tri(any(reshape(sqrt((xy(edges(:,1),1) - xy(edges(:,2),1)).^2 + ...
     (xy(edges(:,1),2) - xy(edges(:,2),2)).^2),[],3) > maxEdge,2),:) = [];

end

% Author:
%       Joseph Kirk
%       jdkirk630@gmail.com
%
function [costs,paths] = Dijkstra(AorV,xyCorE,SID,FID,showWaitbar)
    
    narginchk(2,5);
    
    % Process inputs
    [n,nc] = size(AorV);
    [m,mc] = size(xyCorE);
    if (nargin < 3)
        SID = (1:n);
    elseif isempty(SID)
        SID = (1:n);
    end
    L = length(SID);
    if (nargin < 4)
        FID = (1:n);
    elseif isempty(FID)
        FID = (1:n);
    end
    M = length(FID);
    if (nargin < 5)
        showWaitbar = (n > 1000 && max(L,M) > 1);
    end
    
    
    % Error check inputs
    if (max(SID) > n || min(SID) < 1)
        eval(['help ' mfilename]);
        error('Invalid [SID] input. See help notes above.');
    end
    if (max(FID) > n || min(FID) < 1)
        eval(['help ' mfilename]);
        error('Invalid [FID] input. See help notes above.');
    end
    [E,cost] = process_inputs(AorV,xyCorE);
    
    
    % Reverse the algorithm if it will be more efficient
    isReverse = false;
    if L > M
        E = E(:,[2 1]);
        cost = cost';
        tmp = SID;
        SID = FID;
        FID = tmp;
        isReverse = true;
    end
    
    
    % Initialize output variables
    L = length(SID);
    M = length(FID);
    costs = zeros(L,M);
    paths = num2cell(NaN(L,M));
    
    
    % Create a waitbar if desired
    if showWaitbar
        hWait = waitbar(0,'Calculating Minimum Paths ... ');
    end
    
    
    % Find the minimum costs and paths using Dijkstra's Algorithm
    for k = 1:L
        
        % Initializations
        iTable = NaN(n,1);
        minCost = Inf(n,1);
        isSettled = false(n,1);
        path = num2cell(NaN(n,1));
        I = SID(k);
        minCost(I) = 0;
        iTable(I) = 0;
        isSettled(I) = true;
        path(I) = {I};
        
        % Execute Dijkstra's Algorithm for this vertex
        while any(~isSettled(FID))
            
            % Update the table
            jTable = iTable;
            iTable(I) = NaN;
            nodeIndex = find(E(:,1) == I);
            
            % Calculate the costs to the neighbor nodes and record paths
            for kk = 1:length(nodeIndex)
                J = E(nodeIndex(kk),2);
                if ~isSettled(J)
                    c = cost(I,J);
                    empty = isnan(jTable(J));
                    if empty || (jTable(J) > (jTable(I) + c))
                        iTable(J) = jTable(I) + c;
                        if isReverse
                            path{J} = [J path{I}];
                        else
                            path{J} = [path{I} J];
                        end
                    else
                        iTable(J) = jTable(J);
                    end
                end
            end
            
            % Find values in the table
            K = find(~isnan(iTable));
            if isempty(K)
                break
            else
                % Settle the minimum value in the table
                [~,N] = min(iTable(K));
                I = K(N);
                minCost(I) = iTable(I);
                isSettled(I) = true;
            end
        end
        
        % Store costs and paths
        costs(k,:) = minCost(FID);
        paths(k,:) = path(FID);
        if showWaitbar && ~mod(k,ceil(L/100))
            waitbar(k/L,hWait);
        end
    end
    if showWaitbar
        delete(hWait);
    end
    
    
    % Reformat outputs if algorithm was reversed
    if isReverse
        costs = costs';
        paths = paths';
    end
    
    
    % Pass the path as an array if only one source/destination were given
    if L == 1 && M == 1
        paths = paths{1};
    end
    
    
    % -------------------------------------------------------------------
    function [E,C] = process_inputs(AorV,xyCorE)
        if (n == nc)
            if (m == n)
                A = AorV;
                A = A - diag(diag(A));
                if (m == mc)
                    % Inputs = (A,cost)
                    C = xyCorE;
                    E = a2e(A);
                else
                    % Inputs = (A,xy)
                    xy = xyCorE;
                    E = a2e(A);
                    D = ve2d(xy,E);
                    C = sparse(E(:,1),E(:,2),D,n,n);
                end
            else
                eval(['help ' mfilename]);
                error('Invalid [A,xy] or [A,C] inputs. See help notes above.');
            end
        else
            if (mc == 2)
                % Inputs = (V,E)
                V = AorV;
                E = xyCorE;
                D = ve2d(V,E);
                C = sparse(E(:,1),E(:,2),D,n,n);
            elseif (mc == 3)
                % Inputs = (V,E3)
                E3 = xyCorE;
                E = E3(:,1:2);
                C = sparse(E(:,1),E(:,2),E3(:,3),n,n);
            else
                eval(['help ' mfilename]);
                error('Invalid [V,E] or [V,E3] inputs. See help notes above.');
            end
        end
    end
    
end

% Convert adjacency matrix to edge list
function E = a2e(A)
    [I,J] = find(A);
    E = [I J];
end

% Compute Euclidean distance for edges
function D = ve2d(V,E)
    VI = V(E(:,1),:);
    VJ = V(E(:,2),:);
    D = sqrt(sum((VI - VJ).^2,2));
end


