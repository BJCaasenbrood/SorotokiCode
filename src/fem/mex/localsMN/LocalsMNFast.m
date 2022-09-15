function [Fe,Fb,Qe,Me,Ce,Ke,Kte,Svme,SS,EE,Te,Ve,Vge,Tke,Re,Ue] = ...
        LocalsMNFast(eNode,eDof,dV,Rb,...
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
        MooneyC,...     
        MooneyD...
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
    dJ    = sqrt(det(J0.'*J0));
    
    % deformation gradient   
    F = DeformationGradient(Delta,dNdx,Dim);
    
    % polar decompostion
    [R, Q, ~] = PolarDecomposition(F);

    % increase robustness low density
    Q = Rb*(Q-eye(3)) + eye(3);
    
    % get internal stress matrix
    [S0, D0, Psi] = PiollaStress(MooneyC,MooneyD,F);
    
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

function [S, D, P] = PiollaStress(MooneyC,MooneyD,F)
X12 = 1/2; X13 = 1/3; X23 = 2/3; X43 = 4/3; X53 = 5/3; X89 = 8/9;
MooneyC10 = MooneyC(1);
MooneyC01 = MooneyC(2);
MooneyK   = MooneyD;

C = F.'*F;
C1=C(1,1); C2=C(2,2); C3=C(3,3); C4=C(1,2); C5=C(2,3); C6=C(1,3);
I1 = C1+C2+C3;
I2 = C1*C2+C1*C3+C2*C3-C4^2-C5^2-C6^2;
I3 = det(C);
J1 = I1*I3^(-X13);
J2 = I2*I3^(-X23);
J3 = sqrt(I3);
J3M1 = J3 - 1;
%
I1E = 2*[1,1,1,0,0,0]';
I2E = 2*[C2+C3, C3+C1, C1+C2, -C4, -C5, -C6]';
I3E = 2*[C2*C3-C5^2,  C3*C1-C6^2,  C1*C2-C4^2, ...
    C5*C6-C3*C4, C6*C4-C1*C5, C4*C5-C2*C6]';
%
W1 = I3^(-X13); W2 = X13*I1*I3^(-X43); W3 = I3^(-X23);
W4 = X23*I2*I3^(-X53); W5 = X12*I3^(-X12);
%
J1E = W1*I1E - W2*I3E;
J2E = W3*I2E - W4*I3E;
J3E = W5*I3E;
%
P = MooneyC10*(J1-3) + MooneyC01*(J2-3) + ...
    0.5*MooneyK*(J3 - 1)^2;

Se = MooneyC10*J1E + MooneyC01*J2E + ...
    MooneyK*J3M1*J3E;

S = [Se(1), Se(4), Se(6); 
     Se(4), Se(2), Se(5); 
     Se(6), Se(5), Se(3)];
 
I2EE = [0  4  4  0  0  0; 4  0  4  0  0  0; 4  4  0  0  0  0;
        0  0  0 -2  0  0; 0  0  0  0 -2  0; 0  0  0  0  0 -2];

I3EE = [ 0     4*C3  4*C2  0    -4*C5  0;
         4*C3  0     4*C1  0     0    -4*C6;
         4*C2  4*C1  0    -4*C4  0     0;
         0     0    -4*C4 -2*C3  2*C6  2*C5;
        -4*C5  0     0     2*C6 -2*C1  2*C4;
         0    -4*C6  0     2*C5  2*C4 -2*C2];
%
W1 = X23*I3^(-X12);    W2 = X89*I1*I3^(-X43); W3 = X13*I1*I3^(-X43);
W4 = X43*I3^(-X12);    W5 = X89*I2*I3^(-X53); W6 = I3^(-X23);
W7 = X23*I2*I3^(-X53); W8 = I3^(-X12);        W9 = X12*I3^(-X12);
%
J1EE = -W1*(J1E*J3E' + J3E*J1E') + W2*(J3E*J3E') - W3*I3EE;
J2EE = -W4*(J2E*J3E' + J3E*J2E') + W5*(J3E*J3E') + W6*I2EE - W7*I3EE;
J3EE = -W8*(J3E*J3E') + W9*I3EE;
%
D = MooneyC10*J1EE + MooneyC01*J2EE ...
+ MooneyK*(J3E*J3E') + MooneyK*J3M1*J3EE;
end

%---------------------------------------------------------- polar decompose
function F = DeformationGradient(U,dNdx,Dim)
nn = round(length(U)/Dim);
UU = zeros(nn,Dim);
id1 = round(1:Dim:Dim*nn).';
id2 = round(2:Dim:Dim*nn).';
UU(:,1) = U(id1);
UU(:,2) = U(id2);

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