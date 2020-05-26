classdef YeohMaterial

    properties (Access = public)
        Type = 'Yeoh';
        C1 = 10;
        C2 = 0;
        C3 = 0;
        D1 = 50;
        D2 = 50;
        D3 = 50;
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
function [S, D] = PiollaStress(YeohMaterial,C,R)
%Se = 2nd PK stress [S11, S22, S33, S12, S23, S13];

S = zeros(3,3);
%D = zeros(3,3);
J = sqrt(det(C));

YeohC = [YeohMaterial.C1,YeohMaterial.C2,YeohMaterial.C3];
YeohD = [YeohMaterial.D1,YeohMaterial.D2,YeohMaterial.D3];

I = eye(3,3);
Cinv = minv(C);
C11=C(1,1); 
C22=C(2,2); 
C33=C(3,3);
I1 = C11+C22+C33;
I1iso = J^(-2/3)*I1;
%I1iso = I1;

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

II3 = I - (I1/3)*Cinv;
TOa = TensorOperation(II3,II3,true);
TOb1 = TensorOperation(Cinv,I,true);
TOb2 = TensorOperation(I,Cinv,true);
TOb3 = TensorOperation(Cinv,Cinv,true);
TOb4 = TensorOperation(Cinv,Cinv,false);

TOc = TOb3;
TOd1 = TOb3;
TOd2 = TOb4;

D = (4*J^(-4/3)*alpha)*TOa - ((4/3)*J^(-2/3)*beta)*(TOb1 + TOb2 - ...
    (I1/3)*TOb3 - I1*TOb4) + (J^2)*gamma*TOc + delta*J*(TOd1 -2*TOd2);

end

%---------------------------------------------------------------------- set
function E = Emod(YeohMaterial)
   E = 6*YeohMaterial.C1;
end

end

methods (Access = private)


    
end
end