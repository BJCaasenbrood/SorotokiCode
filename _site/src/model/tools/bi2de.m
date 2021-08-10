function d = bi2de(b, varargin) 
%BI2DE Convert binary vectors to decimal numbers. 
%   D = BI2DE(B) converts a binary vector B to a decimal value D. When B is a 
%   matrix, the conversion is performed row-wise and the output D is a column 
%   vector of decimal values. The default orientation of the binary input 
%   is Right-MSB; the first element in B represents the least significant bit. 
% 
%   In addition to the input matrix, two optional parameters can be given: 
% 
%   D = BI2DE(...,P) converts a base P vector to a decimal value. 
% 
%   D = BI2DE(...,FLAG) uses FLAG to determine the input orientation.  FLAG has 
%   two possible values, 'right-msb' and 'left-msb'.  Giving a 'right-msb' FLAG 
%   does not change the function's default behavior.  Giving a 'left-msb' 
%   FLAG flips the input orientation such that the MSB is on the left. 
% 
%   Examples: 
%   ? B = [0 0 1 1; 1 0 1 0]; 
%   ? T = [0 1 1; 2 1 0]; 
% 
%   ? D = bi2de(B)      ? D = bi2de(B,'left-msb')      ? D = bi2de(T,3) 
%   D =                 D =                            D = 
%       12                   3                             12 
%        5                  10                              5 
% 
%   See also DE2BI. 
 
%   Copyright 1996-2004 The MathWorks, Inc. 
%   $Revision: 1.15.4.1 $  $Date: 2004/08/10 01:32:22 $ 
 
 
% Typical error checking. 
error(nargchk(1,3,nargin)); 
 
% --- Placeholder for the signature string. 
sigStr = ''; 
flag = ''; 
p = []; 
 
% Check the type of the input B 
if ~(isnumeric(b) || islogical(b)) 
    error('The binary input must be numeric or logical.'); 
end 
b = double(b);  % To allow logicals to work 
 
% --- Identify string and numeric arguments 
for i=1:length(varargin) 
   if(i>1) 
      sigStr(size(sigStr,2)+1) = '/'; 
   end 
   % --- Assign the string and numeric flags 
   if(ischar(varargin{i})) 
      sigStr(size(sigStr,2)+1) = 's'; 
   elseif(isnumeric(varargin{i})) 
      sigStr(size(sigStr,2)+1) = 'n'; 
   else 
      error('Optional parameters must be string or numeric.'); 
   end 
end 
 
% --- Identify parameter signitures and assign values to variables 
switch sigStr 
     
    % --- bi2de(d) 
    case '' 
         
	% --- bi2de(d, p) 
	case 'n' 
      p		= varargin{1}; 
 
	% --- bi2de(d, flag) 
	case 's' 
      flag	= varargin{1}; 
 
	% --- bi2de(d, p, flag) 
	case 'n/s' 
      p		= varargin{1}; 
      flag	= varargin{2}; 
 
	% --- bi2de(d, flag, p) 
	case 's/n' 
      flag	= varargin{1}; 
      p		= varargin{2}; 
 
   % --- If the parameter list does not match one of these signatures. 
   otherwise 
      error('Syntax error.'); 
end 
 
if isempty(b) 
   error('Required parameter empty.'); 
end 
 
if max(max(b < 0)) || max(max(~isfinite(b))) || (~isreal(b)) || ... 
     (max(max(floor(b) ~= b))) 
    error('Input must contain only finite real positive integers.'); 
end 
 
% Set up the base to convert from. 
if isempty(p) 
    p = 2; 
elseif max(size(p)) > 1 
   error('Source base must be a scalar.'); 
elseif (floor(p) ~= p) || (~isfinite(p)) || (~isreal(p)) 
   error('Source base must be a finite real integer.'); 
elseif p < 2 
   error('Source base must be greater than or equal to two.'); 
end 
 
if max(max(b)) > (p-1) 
   error('The elements of the matrix are larger than the base can represent.'); 
end 
 
n = size(b,2); 
 
% If a flag is specified to flip the input such that the MSB is to the left. 
if isempty(flag) 
   flag = 'right-msb'; 
elseif ~(strcmp(flag, 'right-msb') || strcmp(flag, 'left-msb')) 
   error('Invalid string flag.'); 
end 
 
if strcmp(flag, 'left-msb') 
 
   b2 = b; 
   b = b2(:,n:-1:1); 
 
end 
 
%%% The conversion 
max_length = 1024; 
pow2vector = p.^(0:1:(size(b,2)-1)); 
size_B = min(max_length,size(b,2)); 
d = b(:,1:size_B)*pow2vector(:,1:size_B).'; 
 
% handle the infs... 
idx = find(max(b(:,max_length+1:size(b,2)).') == 1); 
d(idx) = inf; 
 
% [EOF] 