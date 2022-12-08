clr;
%%
name = 'DragonSkin10';

%%
A = readmatrix([name,'.csv']);

EngineerStrain = A(19:end,4);
EngineerStress = A(19:end,5);
TrueStrain     = A(19:end,2);
TrueStress     = A(19:end,3);

svm = EngineerStress;
lam = EngineerStrain + 1;

svm_ = TrueStress;
lam_ = TrueStrain + 1;

%% plot orignal data
fig(101,[9.25,9]); plot(lam,svm);

%%
[~,lb,ub,Pi0] = YHModel(0,[]);
Model = @(x,Pi) YHModel(x,Pi);
Obj   = @(Pi) Objective(Model,Pi,lam_,svm_);

format short
opt = optimoptions(@fmincon,'FiniteDifferenceStepSize',1e-12);
x   = fmincon(Obj,Pi0,[],[],[],[],lb,ub,[],opt);

disp('Material parameters');
disp(x)

%% AUX functions
% objective
function J = Objective(Model,Pi,Xd,Yd)
[Xe,~] = EnrichStrain(Xd);
ym = Model(Xd,Pi);
ye = Model(Xe,Pi);
Q  = diag(Yd.^-2);
%I  = eye(numel(Yd));
J  = (ym(:) - Yd(:)).'*Q*(ym(:) - Yd(:));

cla;
plot(Xd,Yd,'-','Color',col(10),'LineW',3); hold on;
plot(Xe,ye,'-','Color',col(8)); 
plot(1,1,'Color','w'); 
plot(1,1,'Color','w'); 
plot(1,1,'Color','w'); 

C1 = ['C1 = ', num2str(Pi(1),2)];
C2 = ['C2 = ', num2str(Pi(2),2)];
C3 = ['C3 = ', num2str(Pi(3),2)];

xlim([0, max(Xd)*1.1]);
ylim([-0.2*max(Yd), max(Yd)*1.1]);
legend({'Experiment','Model (Yeoh)',C1,C2,C3},'Location','NorthWest');
xlabel('stretch $\lambda$');
ylabel('engineering stress $\sigma_{11}$');
drawnow;
end

function [y ,y0] = EnrichStrain(x)
y0 = x(:);
y  = [linspace(0.01,1,100).';x(:)];
end

% material models
function [y, lb, ub, Pi0] = YHModel(x,Pi)

    if numel(Pi) == 3
        c1  = Pi(1);
        c2  = Pi(2);
        c3  = Pi(3);
        mat = YeohMaterial([c1,c2,c3]);
    
        y = mat.uniaxial(x);
    else
        y = [];
    end
    
    % bounds and intial
    lb  = [1e-6,-2,-2];
    ub  = [1e6,1e6,1e6];
    Pi0 = [0.,0.0,0.0];
end

function [y, lb, ub, Pi0] = NHModel(x,Pi)

    if numel(Pi) == 2
        E  = Pi(1);
        Nu = Pi(2);
        mat = NeoHookeanMaterial(E,Nu);
    
        y = mat.uniaxial(x);
    else
       y = []; 
    end
    
    % bounds and intial
    lb  = [1e-5,0.1];
    ub  = [25,0.499];
    Pi0 = [0.25,0.3];
end