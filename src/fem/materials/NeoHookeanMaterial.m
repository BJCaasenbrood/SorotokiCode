classdef NeoHookeanMaterial

    properties (Access = public)
        Type = 'NeoHookean';
        E;
        Nu;
        Rho  = 970e-12;
        Zeta = 0.1;
        Lambda;
        Mu; 
    end
    
    properties (Access = private)
        Alpha;
    end

    
%--------------------------------------------------------------------------
methods  
%--------------------------------------------------------------- Mesh Class
function obj = NeoHookeanMaterial(varargin) 
    
    if isempty(varargin)
        E0   = 1;
        Nu0  = 0.45;
    else
        E0 = varargin{1};
        Nu0 = varargin{2};
    end
    
    obj.E      = E0;
    obj.Nu     = Nu0;
    obj.Lambda = (Nu0*E0)/((1+Nu0)*(1-2*Nu0));
    obj.Mu     = E0/(2*(1+Nu0));
end
%---------------------------------------------------------------------- get     
function varargout = get(NeoHookeanMaterial,varargin)
    if nargin > 1
        varargout{nargin-1,1} = [];
        for ii = 1:length(varargin)
            varargout{ii,1} = NeoHookeanMaterial.(varargin{ii});
        end
    else
        varargout = NeoHookeanMaterial.(varargin);
    end
end   
%---------------------------------------------------------------------- get     
function y = getModulus(NeoHookeanMaterial)
y = NeoHookeanMaterial.E;
end
%---------------------------------------------------------------------- set
function NeoHookeanMaterial = set(NeoHookeanMaterial,varargin)
    for ii = 1:2:length(varargin)
        NeoHookeanMaterial.(varargin{ii}) = varargin{ii+1};
    end
end
%------------------------------ 2nd Piolla stress tensor for Neo-Hookean
function [S, D, P] = PiollaStress(NeoHookeanMaterial,F,~)
%Se = 2nd PK stress [S11, S22, S33, S12, S23, S13];

% if (nargin > 2) && Robustness
% %Nu0 = Nu0*0.75;
% end

C100 = NeoHookeanMaterial.Mu/2; 
K    = NeoHookeanMaterial.Lambda/2;

X12 = 1/2; 
X13 = 1/3; 
X23 = 2/3; 
X43 = 4/3; 
X89 = 8/9;
C = F.'*F;

C1=C(1,1); C2=C(2,2); C3=C(3,3); C4=C(1,2); C5=C(2,3); C6=C(1,3);

I1 = C1+C2+C3;
I3 = det(C);
J1 = I1*I3^(-X13);
J3 = sqrt(I3);
J3M1 = J3 - 1;
%
I1E = 2*[1,1,1,0,0,0]';
I3E = 2*[C2*C3-C5^2,  C3*C1-C6^2,  C1*C2-C4^2, ...
    C5*C6-C3*C4, C6*C4-C1*C5, C4*C5-C2*C6]';
%
W1 = I3^(-X13); W2 = X13*I1*I3^(-X43); W5 = X12*I3^(-X12);
%
J1E = W1*I1E - W2*I3E;
J3E = W5*I3E;
%

P = C100*(J1-3) + 0.5*K*(J3 - 1)^2;

Se = C100*J1E + K*J3M1*J3E;

S = [Se(1), Se(4), Se(6); 
     Se(4), Se(2), Se(5); 
     Se(6), Se(5), Se(3)];

I3EE = [ 0     4*C3  4*C2  0    -4*C5  0;
         4*C3  0     4*C1  0     0    -4*C6;
         4*C2  4*C1  0    -4*C4  0     0;
         0     0    -4*C4 -2*C3  2*C6  2*C5;
        -4*C5  0     0     2*C6 -2*C1  2*C4;
         0    -4*C6  0     2*C5  2*C4 -2*C2];
%
W1 = X23*I3^(-X12);    W2 = X89*I1*I3^(-X43); W3 = X13*I1*I3^(-X43);
W8 = I3^(-X12);        W9 = X12*I3^(-X12);
%
J1EE = -W1*(J1E*J3E' + J3E*J1E') + W2*(J3E*J3E') - W3*I3EE;
J3EE = -W8*(J3E*J3E') + W9*I3EE;
%
D = C100*J1EE + K*(J3E*J3E') + K*J3M1*J3EE;

end
%---------------------------------------------------------------------- set
function E = Emod(NeoHookeanMaterial)
   E = NeoHookeanMaterial.E;
end

end
end

