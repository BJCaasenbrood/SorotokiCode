classdef OgdenMaterial

    properties (Access = public)
        Type = 'Ogden';
        C1   = 10;
        C2   = 0;
        A1   = 1;
        A2   = 1;
        Bulk = 1;
        Rho  = 1e-9;
        Zeta = 0.1;
    end
    
    properties (Access = private)
%         ID;
%         SET;
%         WGT;
%         TensorCalc = false;
    end
   
%--------------------------------------------------------------------------
methods  
%--------------------------------------------------------------- Mesh Class
function obj = OgdenMaterial(varargin) 
    
    %[obj.ID, obj.SET, obj.WGT] = Tensor4IdSymmetric;
    
    for ii = 1:2:length(varargin)
        obj.(varargin{ii}) = varargin{ii+1};
    end
end

%---------------------------------------------------------------------- get     
function varargout = get(OgdenMaterial,varargin)
    if nargin > 1
        varargout{nargin-1,1} = [];
        for ii = 1:length(varargin)
            varargout{ii,1} = OgdenMaterial.(varargin{ii});
        end
    else
        varargout = OgdenMaterial.(varargin);
    end
end
        
%---------------------------------------------------------------------- set
function OgdenMaterial = set(OgdenMaterial,varargin)
    for ii = 1:2:length(varargin)
        OgdenMaterial.(varargin{ii}) = varargin{ii+1};
    end
end
    
%------------------------------ 2ND PIOLLA STRESSAND STIFFNESS FOR YEOH
function [S, D, P] = PiollaStress(OgdenMaterial,C)
%Se = 2nd PK stress [S11, S22, S33, S12, S23, S13];
c1 = OgdenMaterial.C1;
c2 = OgdenMaterial.C2;
m1 = OgdenMaterial.A1;
m2 = OgdenMaterial.A2;
kappa = 0*OgdenMaterial.Bulk; % bulk modulus

% spectral decomposition of C
Cd = C(1:2,1:2);%[C(1) C(3); C(3) C(2)];

[Veig,Deig] = eig(Cd);
Deig = diag(Deig);
J2 = det(Cd); % det(C) = J2, where det(F) = J
J = sqrt(J2);

Stretch = Deig.^0.5;
logJ = log(J);

% Energy density W(F) = c1/m1^2*(lambda1^m1+lambda2^m1+lambda3^m1-3-m1*log(J))+c2/m2^2*(lambda1^m2+lambda2^m2+lambda3^m2-3-m2*log(J))+kappa/2*(J-1)^2 (Ogden)
logJ = log(J);
P = c1 / (m1 * m1) * (Stretch(1)^m1 + Stretch(2)^m1 + 1.0 - 3.0 - m1 * logJ) +...
    c2 / (m2 * m2) * (Stretch(1)^m2 + Stretch(2)^m2 + 1.0 - 3.0 - m2 * logJ) +...
    kappa / 2.0 * (J - 1.0) * (J - 1.0);

% Compute derivatives of the energy density with respect to stretches and constants (treating J = lambda1*lambda2*1.0)
dwdL = zeros(2,1);
dwdL(1) = c1 / m1 * Stretch(1)^(m1 - 1.0) + c2 / m2 * Stretch(1)^(m2 - 1.0) - (c1 / m1 + c2 / m2 - J * kappa * (J - 1.0)) / J * Stretch(2);
dwdL(2) = c1 / m1 * Stretch(2)^(m1 - 1.0) + c2 / m2 * Stretch(2)^(m2 - 1.0) - (c1 / m1 + c2 / m2 - J * kappa * (J - 1.0)) / J * Stretch(1);

d2wd2L = zeros(2,2);
d2wd2L(1,1) = c1 / m1 * (m1 - 1.0) * Stretch(1)^(m1 - 2.0) + c2 / m2 * (m2 - 1.0) * Stretch(1)^(m2 - 2.0) + (c1 / m1 + c2 / m2 + J2 * kappa) / J2 * Deig(2);
d2wd2L(2,2) = c1 / m1 * (m1 - 1.0) * Stretch(2)^(m1 - 2.0) + c2 / m2 * (m2 - 1.0) * Stretch(2)^(m2 - 2.0) + (c1 / m1 + c2 / m2 + J2 * kappa) / J2 * Deig(1);
d2wd2L(1,2) = kappa * (2.0 * J - 1.0);
d2wd2L(2,1) = d2wd2L(1,2);

% The second Piola-Kirchhoff stress tensor S
tS = zeros(2);
for i = 1:2
    for j = 1:2
        tS(i,j) = dwdL(1) / Stretch(1) * Veig(i, 1) * Veig(j, 1) +...
            dwdL(2) / Stretch(2) * Veig(i, 2) * Veig(j, 2);
    end
end

% Map tS on Voight notation
idmap = [1, 2, 1; 1, 2, 2]; % map from tensor to Voight notation
S = zeros(3,1);
for i = 1:3
    S(i) = tS(idmap(1,i),idmap(2,i));
end

% Material stiffness D
tD = zeros(2,2,2,2);
for i = 1:2
    for j = 1:2
        for k = 1:2
            for l = 1:2
                tD(i,j,k,l) = 1.0 / Stretch(1) * (-1.0 / Deig(1) * dwdL(1) + 1.0 / Stretch(1) * d2wd2L(1,1)) * Veig(i, 1) * Veig(j, 1) * Veig(k, 1) * Veig(l, 1) +...
                    1.0 / Stretch(1) * (1.0 / Stretch(2) * d2wd2L(1,2)) * Veig(i, 1) * Veig(j, 1) * Veig(k, 2) * Veig(l, 2) +...
                    1.0 / Stretch(2) * (1.0 / Stretch(1) * d2wd2L(2,1)) * Veig(i, 2) * Veig(j, 2) * Veig(k, 1) * Veig(l, 1) +...
                    1.0 / Stretch(2) * (-1.0 / Deig(2) * dwdL(2) + 1.0 / Stretch(2) * d2wd2L(2,2)) * Veig(i, 2) * Veig(j, 2) * Veig(k, 2) * Veig(l, 2);
            end
        end
    end
end
if (abs(Stretch(2) - Stretch(1)) > eps) % two distinct eigenvalues
    for i = 1:2
        for j = 1:2
            for k = 1:2
                for l = 1:2
                    tD(i,j,k,l) = tD(i,j,k,l) + (dwdL(2) / Stretch(2) - dwdL(1) / Stretch(1)) / (Deig(2) - Deig(1)) * (Veig(i, 1) * Veig(j, 2) * Veig(k, 1) * Veig(l, 2) + Veig(i, 1) * Veig(j, 2) * Veig(k, 2) * Veig(l, 1)) +...
                        (dwdL(1) / Stretch(1) - dwdL(2) / Stretch(2)) / (Deig(1) - Deig(2)) * (Veig(i, 2) * Veig(j, 1) * Veig(k, 2) * Veig(l, 1) + Veig(i, 2) * Veig(j, 1) * Veig(k, 1) * Veig(l, 2));
                end
            end
        end
    end
else % multiple eigenvalues
    for i = 1:2
        for j = 1:2
            for k = 1:2
                for l = 1:2
                    tD(i,j,k,l) = tD(i,j,k,l) + 0.5 / Stretch(2) * ((-1.0 / Deig(2) * dwdL(2) + 1.0 / Stretch(2) * d2wd2L(2,2)) - (1.0 / Stretch(1) * d2wd2L(1,2))) * (Veig(i, 1) * Veig(j, 2) * Veig(k, 1) * Veig(l, 2) + Veig(i, 1) * Veig(j, 2) * Veig(k, 2) * Veig(l, 1)) +...
                        0.5 / Stretch(1) * ((-1.0 / Deig(1) * dwdL(1) + 1.0 / Stretch(1) * d2wd2L(1,1)) - (1.0 / Stretch(2) * d2wd2L(2,1))) * (Veig(i, 2) * Veig(j, 1) * Veig(k, 2) * Veig(l, 1) + Veig(i, 2) * Veig(j, 1) * Veig(k, 1) * Veig(l, 2));
                end
            end
        end
    end
end

% Map tD on Voight notation
D = zeros(3,3);
for i = 1:3
    for j = 1:3
        D(i,j) = tD(idmap(1,i),idmap(2,i),idmap(1,j),idmap(2,j));
    end
end

S = [S(1),S(3),0;S(3),S(2),0;0,0,0];
D = [D(1,1),D(1,2),D(1,3),0;D(2,1),D(2,2),0,0;0,0,0,0;0,0,0,D(3,3)];

end

%---------------------------------------------------------------------- set
function E = Emod(YeohMaterial)
   E = 6*YeohMaterial.C1;
end

%------------------------------ 2ND PIOLLA STRESSAND STIFFNESS FOR YEOH
function y = dWdI(YeohMaterial,I1)
C1 = YeohMaterial.C1; 
C2 = YeohMaterial.C2; 
C3 = YeohMaterial.C3;
y = C1 + 2*C2*(I1 -3) + 3*C3*(I1 -3).^2;
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


