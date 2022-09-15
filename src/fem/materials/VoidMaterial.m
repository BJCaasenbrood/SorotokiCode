classdef VoidMaterial

    properties (Access = public)
        Type = 'Void';
        E0   = 0.1;
        C1;
        C2;
        C3;  
        D1   = 1;
        D2   = 1;
        D3   = 1;
        rho  = 1e-2;
        dp   = 0;
        Rho  = 1e-16;
        Zeta = 0.1;
        Rr   = 0.1;
        Cfr  = 1e-6;
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
function obj = VoidMaterial(varargin) 
    
    [obj.ID, obj.SET, obj.WGT] = Tensor4IdSymmetric; 
    
    for ii = 1:2:length(varargin)
        obj.(varargin{ii}) = varargin{ii+1};
    end
    
    obj.C1 = obj.rho*obj.E0;
    obj.C2 = obj.rho*obj.E0*1e-5;
    obj.C3 = obj.rho*obj.E0*1e-5;
    
    obj.D1 = 1e3;
    obj.D2 = 1e3;
    obj.D3 = 1e3;
    
end
%---------------------------------------------------------------------- get     
function varargout = get(VoidMaterial,varargin)
    if nargin > 1
        varargout{nargin-1,1} = [];
        for ii = 1:length(varargin)
            varargout{ii,1} = VoidMaterial.(varargin{ii});
        end
    else
        varargout = VoidMaterial.(varargin);
    end
end       
%---------------------------------------------------------------------- set
function VoidMaterial = set(VoidMaterial,varargin)
    for ii = 1:2:length(varargin)
        VoidMaterial.(varargin{ii}) = varargin{ii+1};
    end
end
%---------------------------------------------------------------------- set
function y = getModulus(VoidMaterial)
    y = 6*VoidMaterial.C1; 
end
%---------------------------------------------------------------------- set
function y = getContactReaction(VoidMaterial)
    y = 6*VoidMaterial.Rr*VoidMaterial.C1; 
end
%---------------------------------------------------------------------- set
function y = getContactFriction(VoidMaterial)
    y = 6*VoidMaterial.Cfr*VoidMaterial.C1; 
end
%------------------------------ 2ND PIOLLA STRESSAND STIFFNESS FOR YEOH
function [S, D, P] = PiollaStress(VoidMaterial,F)

J    = det(F);
Finv = inv(F);   
DP = VoidMaterial.dp;

P0 = zeros(3,3);
for ii = 1:3
    for jj = 1:3
        P0(ii,jj) = DP*J*Finv(jj,ii);
    end
end

S0 = P0*(Finv.');
D00 = zeros(3,3,3,3);

for ii = 1:3
    for jj = 1:3
        for kk = 1:3
            for ll = 1:3
                D00(ii,jj,kk,ll) = DP*J*(Finv(ll,kk)*Finv(jj,ii) - ...
                                 Finv(ll,ii)*Finv(jj,kk));
            end
        end
    end
end

D0 = zeros(6);
Set = [1,1;
       2,2;
       3,3;
       2,3;
       3,1;
       1,2];

for ii = 1:size(Set,1)
    for jj = 1:size(Set,1)
        D0(ii,jj) = D00(Set(ii,1),Set(ii,2),Set(jj,1),Set(jj,2));
    end
end

W = kron([1,2;2,4],ones(3));
D0 = W.*D0;

%-------------------------------------- elastic material
Pe = 0;
Se = zeros(3,3);
C = F.'*F;
J = sqrt(det(C));

VoidC = [VoidMaterial.C1,VoidMaterial.C2,VoidMaterial.C3];
VoidD = [VoidMaterial.D1,VoidMaterial.D2,VoidMaterial.D3];

I = eye(3,3);
Cinv = minv(C);
C11=C(1,1); 
C22=C(2,2); 
C33=C(3,3);
I1 = C11+C22+C33;
I1iso = J^(-2/3)*I1;

for ii = 1:3
    Pe = Pe + VoidC(ii)*(I1iso - 3)^(ii) + (1/VoidD(ii))*(J-1)^(2*ii);
end

for ii = 1:3
    Se = Se + 2*(ii*VoidC(ii)*(I1iso - 3)^(ii-1))*J^(-2/3)*(I - (I1/3)*Cinv)...
        + ((2*ii/VoidD(ii))*(J-1)^(2*ii - 1))*J*Cinv;
end

alpha = 0; beta = 0; gamma = 0; delta = 0;

for ii = 1:3
    kk = ii;
    if ii > 1, alpha = alpha + ii*(ii-1)*VoidC(ii)*(I1iso - 3)^(ii-2); end
    beta = beta + ii*VoidC(ii)*(I1iso - 3)^(ii-1);
    gamma = gamma + (2*kk*(2*kk-1)/VoidD(kk))*(J-1)^(2*kk-2);
    delta = delta + (2*kk/VoidD(kk))*(J-1)^(2*kk-1);
end

II3 = I - (I1/3)*Cinv;
TOa  = TensorOperation(II3,II3,true);
TOb1 = TensorOperation(Cinv,I,true);
TOb2 = TensorOperation(I,Cinv,true);
TOb3 = TensorOperation(Cinv,Cinv,true);
TOb4 = TensorOperation(Cinv,Cinv,false);

TOc = TOb3;
TOd1 = TOb3;
TOd2 = TOb4;

De = (4*J^(-4/3)*alpha)*TOa - ((4/3)*J^(-2/3)*beta)*(TOb1 + TOb2 - ...
    (I1/3)*TOb3 - I1*TOb4) + (J^2)*gamma*TOc + delta*J*(TOd1 -2*TOd2);

% --- combine
D = D0 + De;
P = P0 + Pe;
S = S0 + Se;

end
%---------------------------------------------------------------------- set
function E = Emod(VoidMaterial),  E = 6*VoidMaterial.C1; end


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

function [id, set, W] = Tensor4IdSymmetric

id = [1,1; 2,2; 3,3; 1,2; 2,3; 1,3];

set = {[1,1], [1,2], [1,3], [1,4], [1,5], [1,6],...
       [2,1], [2,2], [2,3], [2,4], [2,5], [2,6],...
       [3,1], [3,2], [3,3], [3,4], [3,5], [3,6],...
       [4,1], [4,2], [4,3], [4,4], [4,5], [4,6],...
       [5,1], [5,2], [5,3], [5,4], [5,5], [5,6],...
       [6,1], [6,2], [6,3], [6,4], [6,5], [6,6]};
   
set = set(:);

%W = kron([1,sqrt(2);sqrt(2),2],ones(3));%%
W = kron([1,2;2,4],ones(3));
% W = [1,1,1,2,2,2;1,1,1,2,2,2;1,1,1,2,2,2;...
%      2,2,2,4,4,4;2,2,2,4,4,4;2,2,2,4,4,4];
end

