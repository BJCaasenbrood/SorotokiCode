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

function X = minv(A)
a=length(A); 
I=eye(a);
augmat=[A I];

for i=1:a-1
    m=augmat(i,i);
    augmat(i,:)=augmat(i,:)/m; %
    for j=i:a-1
        augmat(j+1,:)=augmat(j+1,:) - augmat(i,:)*augmat(j+1,i); 
    end
end
augmat(a,:)=augmat(a,:)/augmat(a,a); 
for k=2:a
    for g=(k-1):-1:1
        augmat(g,:)=augmat(g,:)-augmat(k,:)*augmat(g,k); 
    end
end

X = augmat(:,a+1:2*a); %
end

function T = TensorOperation(A,B,Arg)
%#codegen

a11 = A(1,1); a12 = A(1,2); a13 = A(1,3);  
a21 = A(2,1); a22 = A(2,2); a23 = A(2,3);  
a31 = A(3,1); a32 = A(3,2); a33 = A(3,3);  

b11 = B(1,1); b12 = B(1,2); b13 = B(1,3);  
b21 = B(2,1); b22 = B(2,2); b23 = B(2,3);  
b31 = B(3,1); b32 = B(3,2); b33 = B(3,3);  

if Arg 
    T = [   a11*b11,   a22*b11,   a33*b11, 2*a12*b11, 2*a23*b11, 2*a13*b11;
   a11*b22,   a22*b22,   a33*b22, 2*a12*b22, 2*a23*b22, 2*a13*b22;
   a11*b33,   a22*b33,   a33*b33, 2*a12*b33, 2*a23*b33, 2*a13*b33;
 2*a11*b12, 2*a22*b12, 2*a33*b12, 4*a12*b12, 4*a23*b12, 4*a13*b12;
 2*a11*b23, 2*a22*b23, 2*a33*b23, 4*a12*b23, 4*a23*b23, 4*a13*b23;
 2*a11*b13, 2*a22*b13, 2*a33*b13, 4*a12*b13, 4*a23*b13, 4*a13*b13];
else
    T = [      a11*b11,           a21*b21,           a31*b31,             2*a11*b21,             2*a21*b31,             2*a11*b31;
           a12*b12,           a22*b22,           a32*b32,             2*a12*b22,             2*a22*b32,             2*a12*b32;
           a13*b13,           a23*b23,           a33*b33,             2*a13*b23,             2*a23*b33,             2*a13*b33;
 a11*b12 + a12*b11, a21*b22 + a22*b21, a31*b32 + a32*b31, 2*a11*b22 + 2*a12*b21, 2*a21*b32 + 2*a22*b31, 2*a11*b32 + 2*a12*b31;
 a12*b13 + a13*b12, a22*b23 + a23*b22, a32*b33 + a33*b32, 2*a12*b23 + 2*a13*b22, 2*a22*b33 + 2*a23*b32, 2*a12*b33 + 2*a13*b32;
 a11*b13 + a13*b11, a21*b23 + a23*b21, a31*b33 + a33*b31, 2*a11*b23 + 2*a13*b21, 2*a21*b33 + 2*a23*b31, 2*a11*b33 + 2*a13*b31];
end

end

