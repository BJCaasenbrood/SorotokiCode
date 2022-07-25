function [Fe,Fb,Qe,Me,Ce,Ke,Kte,Svme,SS,EE,Te,Ve,Vge,Tke,Re,Ue] = ...
        LocalsYHFast(eNode,eDof,dV,Rb,...
        Dim,...   % Fem.Dim
        Node0,... % Node0
        Ns,...     % shpfnc for Fem.ShapeFnc{nn}
        dNdxis,... % derive shpfnc for Fem.ShapeFnc{nn}
        W,...     % weights shapefnc
        Utmp,...  % current displacements
        dUtmp,... % current velocities
        Rho,...   % density
        Zeta,...  % dampings
        Grav,...  % gravity
        YeohC,...     
        YeohD...
        )
% get order
nn = length(eNode);
mm = Dim;

% get gauss points and weights
%W   = Fem.ShapeFnc{nn}.W;
Fe  = zeros(mm*nn,1);
Fb  = zeros(mm*nn,1);
Te  = zeros(mm*nn,1);
Me  = zeros(mm*nn,mm*nn);
Ce  = zeros(mm*nn,mm*nn);
Ke  = zeros(mm*nn,mm*nn);
Kte = zeros(mm*nn,mm*nn);
SGP = zeros(length(W),6);
EGP = zeros(length(W),6);
RRe = zeros(3);
UUe = zeros(3);
Qe  = ones(nn,1);
Ve  = 0;
Vge = 0;

Nshp = length(Ns(:,:,1));
NNe  = zeros(length(W),Nshp);

if mm == 2, Et = [dV;dV;0];
else, Et = [dV;dV;dV;0;0;0];
end

% get displacement field
Delta  = Utmp(eDof,:);
dDelta = dUtmp(eDof,:);

% quadrature loop
for q = 1:length(W)
    
    % extract shape-functions
    N     = Ns(:,:,q);
    dNdxi = dNdxis(:,:,q);
    J0    = Node0(eNode,:).'*dNdxi;
    dNdx  = dNdxi/J0;
    dJ    = det(J0);
    
    % deformation gradient   
    F = DeformationGradient(Delta,dNdx,Dim);
    
    % polar decompostion
    [R, Q, ~] = PolarDecomposition(F);

    % increase robustness low density
    Q = Rb*(Q-eye(3)) + eye(3);
    
    % get internal stress matrix
    [S0, D0, Psi] = PiollaStress(YeohC,YeohD,F);
    
    % voigt-notation vectorize
    S = VoightNotation(S0);
    
    % reduced isotropic matrices
    [Se, De, Ge] = IsotropicReduction(D0,S,Dim);
    
    % nonlinear strain-displacement operator
    [Bnl,Bg,NN,tau] = NonlinearStrainOperatorFast(N,dNdx,F);
    
    % local elemental rotation
    RRe = RRe + R/nn;
    UUe = UUe + Q/nn;
    
    % internal force vector
    Fe = Fe + tau*W(q)*Bnl.'*Se*dJ;
    
    % (graviational) body force vector
    Fb = Fb + tau*W(q)*Rho*(NN.')*Grav(:)*dJ;
    
    % lineararized stiffness matrix
    Ke = Ke + tau*W(q)*(Bnl.'*De*Bnl)*dJ;
    
    % tangent stiffness matrix
    Kte = Kte + tau*W(q)*(Bnl.'*De*Bnl + Bg.'*Ge*Bg)*dJ;
    
    % mass matrix
    Me = Me + tau*W(q)*Rho*(NN.')*NN*dJ;
    
    % dampings matrix
    Ce = Ce + Zeta*Me;
    
    % thermal expansion force
    Te = Te + tau*W(q)*Bnl.'*De*Et*dJ;
    
    % elemental potential energy
    Ve = Ve  + tau*W(q)*Psi*dJ;
    
    % graviational energy
    Vge = Vge - tau*W(q)*((N.'*Node0(eNode,:)) + ...
        (NN*Delta).')*Rho*Grav(:)*dJ;
    
    % lagrangian strain
    Elagran = (1/2)*(F.'*F - eye(3));
    
    % true stress and lagrangian strain
    SGP(q,:) = VoightNotation((1/det(F))*F*S0*(F.'));
    EGP(q,:) = VoightNotation(Elagran);
    
    % construct shape functions
    NNe(((q-1)*Nshp + 1):(q*Nshp)) = N(:).';
end

% compute elemental kinetic energy
Tke = 0.5*(dDelta).'*Me*(dDelta);    

%compute elemental rotation matrix
[Ur,~,Vr] = svd(RRe);
Re = (Ur*Vr.');
Ue = UUe;

SS = NNe.'*SGP;
EE = NNe.'*EGP;
[Svm, ~] = VonMises(SS(:,1),SS(:,2),SS(:,3),...
                    SS(:,4),SS(:,5),SS(:,6));
Svme = Svm(:); 

end

function [S, D, P] = PiollaStress(YeohC,YeohD,F)
%Se = 2nd PK stress [S11, S22, S33, S12, S23, S13];
P = 0;
S = zeros(3,3);
C = F.'*F;
J = sqrt(det(C));

% YeohC = [YeohMaterial.C1,YeohMaterial.C2,YeohMaterial.C3];
% YeohD = [YeohMaterial.D1,YeohMaterial.D2,YeohMaterial.D3];

I = eye(3,3);
Cinv = minv(C);
C11=C(1,1); 
C22=C(2,2); 
C33=C(3,3);
I1 = C11+C22+C33;
I1iso = J^(-2/3)*I1;
%I1iso = I1;

for ii = 1:3
    P = P + YeohC(ii)*(I1iso - 3)^(ii) + (1/YeohD(ii))*(J-1)^(2*ii);
end

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
TOa  = TensorOperation(II3,II3,true);
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

%------------------------------ 2ND PIOLLA STRESSAND STIFFNESS FOR YEOH
function y = dWdI(YeohMaterial,I1)
c1 = YeohMaterial.C1; 
c2 = YeohMaterial.C2; 
c3 = YeohMaterial.C3;
y = c1 + 2*c2*(I1 -3) + 3*c3*(I1 -3).^2;
end

%------------------------------ matrix inverse
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

%------------------------------ tensor operations
function T = TensorOperation(A,B,Arg)

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

%---------------------------------------------------------- polar decompose
function F = DeformationGradient(U,dNdx,Dim)
nn = length(U)/Dim;
UU = zeros(nn,Dim);
UU(:,1) = U(1:Dim:Dim*nn);
UU(:,2) = U(2:Dim:Dim*nn);

if Dim == 2
    F = (dNdx'*UU)';
    F = [F(1,1)+1,F(1,2),0; F(2,1),F(2,2)+1,0;0,0,1];
else
    UU(:,3) = U(3:Dim:Dim*nn);
    F = (dNdx'*UU)' + eye(3);
end
end
%---------------------------------------------------------- polar decompose
function [R,S,V] = PolarDecomposition(F)
C = F.'*F;
[Q0, lambdasquare] = eig(C);

lambda = sqrt(diag((lambdasquare))); 
Uinv   = repmat(1./lambda',size(F,1),1).*Q0*Q0.';

R = real(F*Uinv);
S = real(R.'*F);
V = real(F*R.');
end
%------------------------------------------------ nonlinear strain operator
function [Bn,Bg,NN,tau] = NonlinearStrainOperatorFast(N,dNdx,F)
nn = length(N);
mm = size(dNdx,2);
zz = mm*nn;

NN = zeros(mm,zz);

id1 = 1:mm:zz;
id2 = 2:mm:zz;

NN(1,id1) = N.';
NN(2,id2) = N.';
dNdxX = dNdx(:,1).';
dNdxY = dNdx(:,2).';

if mm == 3
    NN(3,3:mm:zz) = N.';
    dNdxZ = dNdx(:,3).';
else
    dNdxZ = 0;
end

Bn = zeros((mm-1)*3,zz);
Bg = zeros((mm-1)*4+(mm-2),zz);

if mm == 2 % 2-dimensional
    Bn(1,id1) = dNdxX*F(1,1);
    Bn(1,id2) = dNdxX*F(2,1);
    Bn(2,id1) = dNdxY*F(1,2);
    Bn(2,id2) = dNdxY*F(2,2);
    Bn(3,id1) = dNdxX*F(1,2) + dNdxY*F(1,1);
    Bn(3,id2) = dNdxX*F(2,2) + dNdxY*F(2,1);
    
    Bg(1,id1) = dNdxX;
    Bg(2,id1) = dNdxY;
    Bg(3,id2) = dNdxX;
    Bg(4,id2) = dNdxY;
    
else % 3-dimensional
    Bn(1,1:mm:mm*nn) = dNdxX*F(1,1);
    Bn(1,2:mm:mm*nn) = dNdxX*F(2,1);
    Bn(1,3:mm:mm*nn) = dNdxX*F(3,1);
    Bn(2,1:mm:mm*nn) = dNdxY*F(1,2);
    Bn(2,2:mm:mm*nn) = dNdxY*F(2,2);
    Bn(2,3:mm:mm*nn) = dNdxY*F(3,2);
    Bn(3,1:mm:mm*nn) = dNdxZ*F(1,3);
    Bn(3,2:mm:mm*nn) = dNdxZ*F(2,3);
    Bn(3,3:mm:mm*nn) = dNdxZ*F(3,3);
    Bn(4,1:mm:mm*nn) = dNdxX*F(1,2) + dNdxY*F(1,1);
    Bn(4,2:mm:mm*nn) = dNdxX*F(2,2) + dNdxY*F(2,1);
    Bn(4,3:mm:mm*nn) = dNdxX*F(3,2) + dNdxY*F(3,1);
    Bn(5,1:mm:mm*nn) = dNdxY*F(1,3) + dNdxZ*F(1,2);
    Bn(5,2:mm:mm*nn) = dNdxY*F(2,3) + dNdxZ*F(2,2);
    Bn(5,3:mm:mm*nn) = dNdxY*F(3,3) + dNdxZ*F(3,2);
    Bn(6,1:mm:mm*nn) = dNdxX*F(1,3) + dNdxZ*F(1,1);
    Bn(6,2:mm:mm*nn) = dNdxX*F(2,3) + dNdxZ*F(2,1);
    Bn(6,3:mm:mm*nn) = dNdxX*F(3,3) + dNdxZ*F(3,1);
    
    Bg(1,1:mm:mm*nn) = dNdxX;
    Bg(2,1:mm:mm*nn) = dNdxY;
    Bg(3,1:mm:mm*nn) = dNdxZ;
    Bg(4,2:mm:mm*nn) = dNdxX;
    Bg(5,2:mm:mm*nn) = dNdxY;    
    Bg(6,2:mm:mm*nn) = dNdxZ;   
    Bg(7,3:mm:mm*nn) = dNdxX;
    Bg(8,3:mm:mm*nn) = dNdxY;    
    Bg(9,3:mm:mm*nn) = dNdxZ;     
end

tau = 1;

end 
%------------------------------------------------ nonlinear strain operator
function [S, D, G] = IsotropicReduction(D0,S0,Dim)

if Dim == 2
    G = zeros(4,4);
    D   = [D0(1,1), D0(1,2),       0;
           D0(2,1), D0(2,2),       0;
                 0,       0, D0(4,4)];
    SIG = [S0(1), S0(4); S0(4), S0(2)];
    S   = [S0(1); S0(2); S0(4)]; 
    
    G(1:2,1:2) = SIG;
    G(3:4,3:4) = SIG;
else
    G = zeros(9,9);
    D = D0;
    S = S0;
    SIG = [S0(1), S0(4), S0(6);
           S0(4), S0(2), S0(5);
           S0(6), S0(5),S0(3)];
       
    G(1:3,1:3) = SIG;
    G(4:6,4:6) = SIG;
    G(7:9,7:9) = SIG;
end

end
%------------------------------------------------ nonlinear strain operator
function [Svm, Svmm] = VonMises(S11,S22,S33,S12,S23,S13)
s11 = S11; s22 = S22; s33 = S33; s12 = S12; s23 = S23; s13 = S13;
Svm = sqrt(0.5*((s11-s22).^2 + (s22-s33).^2 + (s33-s11).^2 ...
    + 6*(s12.^2 + s23.^2 + s13.^2)));
Svmm = mean(Svm); 
end
%--------------------------------------------------- Kelvin-voight notation
function Sv = VoightNotation(S)
Sv = [S(1,1); S(2,2); S(3,3); S(1,2); S(2,3); S(1,3)]; 
end