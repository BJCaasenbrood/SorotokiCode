classdef YeohMaterial

    properties (Access = public)
        Type = 'Yeoh';
        C1 = 10;
        C2 = 0;
        C3 = 0;
        D1 = 500;
        D2 = 500;
        D3 = 500;
    end
    
    properties (Access = private)
        ID;
        SET;
        WGT;
        TensorCalc = false;
    end
   
%--------------------------------------------------------------------------
methods  
%--------------------------------------------------------------- Mesh Class
function obj = YeohMaterial(varargin) 
    
    [obj.ID, obj.SET, obj.WGT] = Tensor4IdSymmetric;
    
    for ii = 1:2:length(varargin)
        obj.(varargin{ii}) = varargin{ii+1};
    end
end

%---------------------------------------------------------------------- get     
function varargout = get(YeohMaterial,varargin)
    if nargin > 1
        varargout{nargin-1,1} = [];
        for ii = 1:length(varargin)
            varargout{ii,1} = YeohMaterial.(varargin{ii});
        end
    else
        varargout = YeohMaterial.(varargin);
    end
end
        
%---------------------------------------------------------------------- set
function YeohMaterial = set(YeohMaterial,varargin)
    for ii = 1:2:length(varargin)
        YeohMaterial.(varargin{ii}) = varargin{ii+1};
    end
end
    
%------------------------------ 2ND PIOLLA STRESSAND STIFFNESS FOR YEOH
function [S, D] = PiollaStress(YeohMaterial,C,Robustness)
%Se = 2nd PK stress [S11, S22, S33, S12, S23, S13];

S = zeros(3,3);
%D = zeros(3,3);
J = sqrt(det(C));
YeohC = [YeohMaterial.C1,YeohMaterial.C2,YeohMaterial.C3];
YeohD = [YeohMaterial.D1,YeohMaterial.D2,YeohMaterial.D3];
%Fvol = J^(1/3)*eye(3);
%Fiso = J^(-1/3)*F;

if nargin > 2
YeohC(2) = YeohC(2)*Robustness^3;
end

I = eye(3,3);
Cinv = C\I;
C11=C(1,1); C22=C(2,2); C33=C(3,3);
I1 = C11+C22+C33;
I1iso = J^(-2/3)*I1;

for ii = 1:3
    S = S + 2*(ii*YeohC(ii)*(I1iso - 3)^(ii-1))*J^(-2/3)*(I - (I1/3)*Cinv)...
        + ((2*ii/YeohD(ii))*(J-1)^(2*ii - 1))*J*Cinv;
end

alpha = 0; beta = 0; gamma = 0; delta = 0;

for ii = 1:3
    kk = ii;
    if ii > 1, alpha = alpha + ii*(ii-1)*YeohC(ii)*(I1iso - 3)^(ii-2); end
    beta = beta + ii*YeohC(ii)*(I1iso - 3)^(ii-1);
    gamma = gamma + (2*kk*(2*kk-1)/YeohD(kk))*(J-1)^(2*kk-2);
    delta = delta + (2*kk/YeohD(kk))*(J-1)^(2*kk-1);
end

% TOa = TensorOperation(YeohMaterial,I - (I1/3)*Cinv,I - (I1/3)*Cinv,'x');
% TOb1 = TensorOperation(YeohMaterial,Cinv,I,'x');
% TOb2 = TensorOperation(YeohMaterial,I,Cinv,'x');
% TOb3 = TensorOperation(YeohMaterial,Cinv,Cinv,'x');
% TOb4 = TensorOperation(YeohMaterial,Cinv,Cinv,'xt');
TOa = TensorOperation(I - (I1/3)*Cinv,I - (I1/3)*Cinv,'x');
TOb1 = TensorOperation(Cinv,I,'x');
TOb2 = TensorOperation(I,Cinv,'x');
TOb3 = TensorOperation(Cinv,Cinv,'x');
TOb4 = TensorOperation(Cinv,Cinv,'xt');

TOc = TOb3;
TOd1 = TOb3;
TOd2 = TOb4;

D = (4*J^(-4/3)*alpha)*TOa + ((-4/3)*J^(-2/3)*beta)*(TOb1 + TOb2 - ...
    (I1/3)*TOb3 - I1*TOb4) + (J^2)*gamma*TOc + delta*J*(TOd1 -2*TOd2);

end

end

methods (Access = private)

% %-------------------------------------- TENSOR PERMUTATION SETS FOR 3:3/3:3
% function T = TensorOperation(YeohMat,A,B,Arg)
% 
% id = YeohMat.ID; 
% set = YeohMat.SET; 
% W = YeohMat.WGT; 
% 
% T = zeros(length(id)^2,1);
% 
% if strcmp(Arg,'x') % kronecker product
%     for kk = 1:length(set)
%         row = set{kk}(1);
%         col = set{kk}(2);
%         i = id(row,1);
%         j = id(row,2);
%         k = id(col,1);
%         l = id(col,2);
%         Aij = A(i,j);
%         Bkl = B(k,l);
%         T(kk) = Aij*Bkl;
%     end
% elseif strcmp(Arg,'xt') % symmetric kronecker product
%     for kk = 1:length(set)
%         row = set{kk}(1);
%         col = set{kk}(2);
%         i = id(row,1);
%         j = id(row,2);
%         k = id(col,1);
%         l = id(col,2);
%         Aik = A(i,k);
%         Ail = A(i,l);
%         Bjk = B(j,k);
%         Bjl = B(j,l);
%         
%         T(kk) = 0.5*(Aik*Bjl + Ail*Bjk);
%     end
% end
% 
% T = reshape(T,length(id),length(id));
% T = T.*W;
% end

end
end

