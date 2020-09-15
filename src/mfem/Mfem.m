function x = MFAnalysis(x,Request,Arg)
if strcmp(x,'Help') || strcmp(x,'help')
    Help(); return;
else
switch(Request)
    case('Inductance'); x = computeInductanceNetwork(x);
    case('Bfield');     x = generateBField(x,Arg);
end
end

end

%--------------------------------------------------------------------- HELP
function x = Help()
x = '';
CallDisplay('');
CallDisplay('-----------------------------------------------------------');
CallDisplay('function: x = MFAnaylsis(x,Request,Arg)');
CallDisplay('-----------------------------------------------------------');

info = ...
    [' BLENDER is system of functions to transform or deform FV',...
    '\n \t vertice-based meshes in R3'];

CallDisplay(info);
end

%--------------------------------------------------------------------- HELP
function mag = computeInductanceNetwork(mag)

NElem = mag.NElem;
Element = mag.Element;
Node = mag.Node;
mu0 = mag.mu;
Y = mag.Y;

dcf = NElem*1e3;

if ~isfield(mag,'dC')
[dC,dR,Cm] = LocalElement(Node,Element);
mag.dC = dC;
mag.dR = dR;
mag.Cm = Cm;
else
dC = mag.dC;
dR = mag.dR;
Cm = mag.Cm;  
end

PermMat = PermutationMatrix(NElem);
if ~PermMat, return; end

CallDisplay('performing neumann integral...')
CallExecuted(['number of combination = ',num2str(0.5*NElem*(NElem-1))]);

Ltmp = zeros(NElem,NElem); Iter = 1; Rmin = 1e-6;%min(dR);

for k = 1:length(PermMat)
    id = PermMat(k,1); jd = PermMat(k,2);
    
    Rtmp = norm(Cm(id,:)-Cm(jd,:));
    
    % check if 
    if Rtmp > Rmin/2
        Ltmp(id,jd) = dot(dC(id,:),dC(jd,:).')/Rtmp;
        Ltmp(jd,id) = Ltmp(id,jd);
    end
    
    if ~mod(Iter,dcf)
        CallExecuted(['batch: ',num2str((Iter-dcf)/1e3),'k','-',...
            num2str(Iter/1e3),'k']);
    end
    
    Iter = Iter + 1;
end

mag.length = sum(dR);
mag.inductance = 0.25*mu0/pi*((sum(sum(Ltmp)) + sum(dR)*Y));
mag.inductance_ = mu0*sum(sum(Ltmp))/(pi*4);
end

%--------------------------------------------------------------------- HELP
function mag = generateBField(mag,N)
%if nargin == 1, N = 10; end
mu0 = 4*pi*1e-7; 
BdBox = mag.BdBox;

x = linspace(BdBox(1),BdBox(2),N(1));
y = linspace(BdBox(3),BdBox(4),N(2)); 
z = linspace(BdBox(5),BdBox(6),N(3));

[X,Y,Z] = meshgrid(x,y,z);

BX = zeros(size(X,1),size(X,2),size(X,3));
BY = zeros(size(X,1),size(X,2),size(X,3));
BZ = zeros(size(X,1),size(X,2),size(X,3));

for nF = 1%:mag.N % Loop on each coil

    C = mag.Node;
    %dC = generateWireLengths(C);
    I = 1;
    
    % Discretization of C
    x_P = C(:,1).'; 
    y_P = C(:,2).'; 
    z_P = C(:,3).'; 
    
    DX = size(C,1)-1; % Number of points defining C
    
    for i = 1:DX % Loop on the segments defining C
        L_C_i = norm(C(i,:)-C(i+1,:));
        NP = 10; % Number of points required to have a discretization step smaller than dC
        x_P = [x_P,linspace(C(i,1), C(i+1,1), NP)]; % discretization of C for x component
        y_P = [y_P,linspace(C(i,2), C(i+1,2), NP)]; % discretization of C for y component
        z_P = [z_P,linspace(C(i,3), C(i+1,3), NP)]; % discretization of C for z component
    end

    % add contribution of each source point P on each field point M (where we want to calculate the field)
    for m = 1:size(X,1)
        for n = 1:size(X,2)
            for p = 1:size(X,3)

            % M is the field point
            x_M = X(m,n,p);
            y_M = Y(m,n,p);
            z_M = Z(m,n,p);
            %DBx = 0; DBy = 0; DBz = 0;
            % loop on each discretized segment of C PkPk+1
            for k = 1:length(x_P)-1
                PkM3 = (sqrt((x_M-x_P(k))^2 + (y_M-y_P(k))^2 + (z_M-z_P(k))^2))^3;
                DBx(k) = ((y_P(k+1)-y_P(k))*(z_M-z_P(k))-(z_P(k+1)-z_P(k))*(y_M-y_P(k)))/PkM3;
                DBy(k) = ((z_P(k+1)-z_P(k))*(x_M-x_P(k))-(x_P(k+1)-x_P(k))*(z_M-z_P(k)))/PkM3;
                DBz(k) = ((x_P(k+1)-x_P(k))*(y_M-y_P(k))-(y_P(k+1)-y_P(k))*(x_M-x_P(k)))/PkM3;
            end
            
            % summation
            BX(m,n,p) = BX(m,n,p) + mu0*I/4/pi*sum(DBx);
            BY(m,n,p) = BY(m,n,p) + mu0*I/4/pi*sum(DBy);
            BZ(m,n,p) = BZ(m,n,p) + mu0*I/4/pi*sum(DBz);

            end
        end
    end
end

mag.X = X;
mag.Y = Y;
mag.Z = Z;
mag.BX = BX;
mag.BY = BY;
mag.BZ = BZ;

end

%--------------------------------------------------------------------------
function [dC,dR,Cm] = LocalElement(Node,Element)
SubElem = num2cell(Element,2);
dC = cellfun(@(E) diff(Node(E,:)),SubElem,'UniformOutput',false);
Cm = cellfun(@(E) mean(Node(E,:)),SubElem,'UniformOutput',false);
dR = cellfun(@(E) norm(E,2),dC,'UniformOutput',false);
%find the relative segment components

dC = vertcat(dC{:}); 
Cm = vertcat(Cm{:}); 
dR = vertcat(dR{:}); 
end

%--------------------------------------------------------------------------
function PMat = PermutationMatrix(NElem)
usermem = memory; MaxArray = usermem.MaxPossibleArrayBytes / 8;
if (NElem*NElem > MaxArray )
    CallWarning('not enough dynamic memory available!','MFAnalysis');
    CallWarning('reduce the number of elements','MFAnalysis');
    PMat = 0;    
else
    %PMat = sparse(0.5*NElem*(NElem-1),2); 
    PMat = VChooseK(1:NElem,2);
end
end