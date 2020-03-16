classdef NeoHookeanMaterial

    properties (Access = public)
        Type = 'NeoHookean';
        E = 10;
        Nu = 1;
        C10;
        D1 = 10;
    end
    
    properties (Access = private)
        Lambda
        Mu; 
    end

    
%--------------------------------------------------------------------------
methods  
%--------------------------------------------------------------- Mesh Class
function obj = NeoHookeanMaterial(varargin) 
    
    for ii = 1:2:length(varargin)
        obj.(varargin{ii}) = varargin{ii+1};
    end
    
    if isempty(obj.C10)
    E0 = obj.E;
    Nu0 = obj.Nu;
    obj.Lambda = (Nu0*E0)/((1+Nu0)*(1-2*Nu0));
    obj.Mu = E0/(2*(1+Nu0));
    obj.C10 = obj.Mu/2;
    obj.D1 = obj.Lambda/2;
    end
    
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
        
%---------------------------------------------------------------------- set
function NeoHookeanMaterial = set(NeoHookeanMaterial,varargin)
    for ii = 1:2:length(varargin)
        NeoHookeanMaterial.(varargin{ii}) = varargin{ii+1};
    end
end

%------------------------------ 2ND PIOLLA STRESSAND STIFFNESS FOR YEOH
function [S, D] = PiollaStress(NeoHookeanMaterial,C,Robustness)
%Se = 2nd PK stress [S11, S22, S33, S12, S23, S13];

if (nargin > 2) && Robustness
%Nu0 = Nu0*0.75;
end

C100 = NeoHookeanMaterial.C10; K = NeoHookeanMaterial.D1;

X12 = 1/2; X13 = 1/3; X23 = 2/3; X43 = 4/3; X89 = 8/9;

C1=C(1,1); C2=C(2,2); C3=C(3,3); C4=C(1,2); C5=C(2,3); C6=C(1,3);
I1 = C1+C2+C3;
I2 = C1*C2+C1*C3+C2*C3-C4^2-C5^2-C6^2;
I3 = det(C);
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

