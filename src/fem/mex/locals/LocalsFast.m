function [Fe,Fb,Qe,Me,Ce,Ke,Kte,Svme,SS,EE,Te,Ve,Vge,Tke,Re,Ue] = ...
        LocalsFast(Fem,eNode,eDof,dV,Rb)
% get order
nn = length(eNode);
mm = Fem.Dim;

% get gauss points and weights
W   = Fem.ShapeFnc{nn}.W;
%Q   = Fem.ShapeFnc{nn}.Q;    
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
%Tke = 0;

Nshp = length(Fem.ShapeFnc{nn}.N(:,:,1));
NNe  = zeros(length(W),Nshp);

if mm == 2
    Et = [dV;dV;0];
else
    Et = [dV;dV;dV;0;0;0];
end

% get displacement field
Delta  = Fem.Utmp(eDof,:);
dDelta = Fem.dUtmp(eDof,:);

% computation of Piolla stress
PiollaStress = @(x) Fem.Material.PiollaStress(x);

% quadrature loop
for q = 1:length(W)
    
    % extract shape-functions
    N     = Fem.ShapeFnc{nn}.N(:,:,q);
    dNdxi = Fem.ShapeFnc{nn}.dNdxi(:,:,q);
    J0    = Fem.Node0(eNode,:).'*dNdxi;
    dNdx  = dNdxi/J0;
    dJ    = det(J0);
    
    % deformation gradient   
    F = DeformationGradient(Fem,Delta,dNdx);
    
    % polar decompostion
    [R, Q, ~] = PolarDecomposition(Fem,F);
    %[R, Q, ~] = PolarDecompositionFast(F);

    % increase robustness low density
    Q = Rb*(Q-eye(3)) + eye(3);
    
    % reconstruct deformation gradient
    %F = R*Q;
    
%     % right cauchy-green strain
%     C = F.'*F;
%     
%     if strcmp(Fem.Type,'PlaneStress') && Fem.Dim < 3
%         C(3,3) = det(F)/(C(1,1)*C(2,2) - C(1,2)*C(2,1));
%     end

    % get internal stress matrix
    [S0, D0, Psi] = PiollaStress(F);
    
    % voigt-notation vectorize
    S = VoightNotation(S0);
    
    % reduced isotropic matrices
    [Se, De, Ge] = IsotropicReduction(Fem,D0,S);
    
    % nonlinear strain-displacement operator
    %[Bnl,Bg,NN,tau] = NonlinearStrainOperator(Fem,N,dNdx,F);
    [Bnl,Bg,NN,tau] = NonlinearStrainOperatorFast(N,dNdx,F);
    
    % local elemental rotation
    RRe = RRe + R/nn;
    UUe = UUe + Q/nn;
    
    % internal force vector
    Fe = Fe + tau*W(q)*Bnl.'*Se*dJ;
    
    % (graviational) body force vector
    Fb = Fb + tau*W(q)*Fem.Material.Rho*(NN.')*Fem.Gravity(:)*dJ;
    
    % lineararized stiffness matrix
    Ke = Ke + tau*W(q)*(Bnl.'*De*Bnl)*dJ;
    
    % tangent stiffness matrix
    Kte = Kte + tau*W(q)*(Bnl.'*De*Bnl + Bg.'*Ge*Bg)*dJ;
    
    % mass matrix
    Me = Me + tau*W(q)*Fem.Material.Rho*(NN.')*NN*dJ;
    
    % dampings matrix
    Ce = Ce + Fem.Material.Zeta*Me;
    
    % thermal expansion force
    Te = Te + tau*W(q)*Bnl.'*De*Et*dJ;
    
    % elemental potential energy
    Ve = Ve  + tau*W(q)*Psi*dJ;
    
    % graviational energy
    Vge = Vge - tau*W(q)*((N.'*Fem.Node0(eNode,:)) + ...
        (NN*Delta).')*Fem.Material.Rho*Fem.Gravity(:)*dJ;
    
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