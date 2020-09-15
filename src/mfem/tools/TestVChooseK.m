function TestVChooseK(doSpeed)
% Automatic test: VChooseK
% This is a routine for automatic testing. It is not needed for processing and
% can be deleted or moved to a folder, where it does not bother.
%
% TestVChooseK(doSpeed)
% INPUT:
%   doSpeed: Optional logical flag to trigger time consuming speed tests.
%            Default: TRUE. If no speed test is defined, this is ignored.
% OUTPUT:
%   On failure the test stops with an error.
%
% The test can be really slow, if less than 750 MB free RAM is available.
% Either drink a cup of coffee, buy more RAM or start this function with the
% argument 0 to reduce the speed testing.
%
% Some other programs from the FEX are called, if they are found. Some are not
% compatible with Matlab 6.5 and "crash" is displayed for the speed
% measurements.
%
% Tested: Matlab 6.5, 7.7, 7.8, WinXP
% Author: Jan Simon, Heidelberg, (C) 2009 matlab.THIS_YEAR(a)nMINUSsimonDOTde
% License: BSD (use/copy/modify on own risk, but mention author)

% $JRev: R0i V:007 Sum:GknnodeJlhUz Date:24-Dec-2008 02:59:40 $
% $File: User\JSim\Published\VChooseK\TestVChooseK.m $

% Initialize: ==================================================================
ErrID = ['JSim:', mfilename];
LF = char(10);

% Default for the input:
if nargin == 0
   doSpeed = true;
end

% Include the local M version VChooseK_M in the speed test:
% (actually not needed, but might be interesting)
testLocalMVersion = true;

% Times for testing:
if doSpeed
   randDelay = 5.0;    % [sec], time for random tests
   SpeedTime = 1.0;    % [sec], time for speed tests
else
   randDelay = 2.0;
   SpeedTime = 0.25;
end

% Hello:
whichFunc = which('VChooseK');
disp(['==== Test VChooseK:  ', datestr(now, 0), LF, ...
      'Version: ', whichFunc, LF]);

% Start tests: -----------------------------------------------------------------
% Loop over accepted classes:
KList     = [1, 2, 3, 5, 7, 10];
ClassList = {'double', 'single', 'uint32', 'int16', 'int8', 'char'};
for K = KList
   for iClass = 1:length(ClassList)
      aClass = ClassList{iClass};
      disp(['== ', upper(aClass), ', ', sprintf('K = %d:', K)]);
      
      % Standard tests - empty input, cell without populated elements, small
      % known answer test:
      S = VChooseK(my_cast([], aClass), K);
      if isempty([]) && isa(S, aClass)
         disp('  ok: empty matrix');
      else
         error(ErrID, 'Failed for empty matrix.');
      end
      
      S = VChooseK(my_cast(1:K, aClass), K);
      if isequal(S, 1:K) && isa(S, aClass)
         disp('  ok: 1:K');
      else
         error(ErrID, 'Failed for 1:K.');
      end
      
      S = VChooseK(my_cast(1:K+1, aClass), K);
      T = nchoosek(my_cast(1:K+1, aClass), K);
      if isequal(S, T) && isa(S, aClass)
         disp('  ok: 1:K+1');
      else
         error(ErrID, 'Failed for 1:K+1.');
      end
      
      S = VChooseK(my_cast(2:K+3, aClass), K);
      T = nchoosek(my_cast(2:K+3, aClass), K);
      if isequal(S, T) && isa(S, aClass)
         disp('  ok: 2:K+3');
      else
         error(ErrID, 'Failed for 2:K+3.');
      end
      
      S = VChooseK(my_cast(transpose(3:K+4), aClass), K);
      T = VChooseK(my_cast(3:K+4, aClass), K);
      if isequal(S, T) && isa(S, aClass)
         disp('  ok: transpose(3:K+4)');
      else
         error(ErrID, 'Failed for transpose(3:K+4).');
      end
   end  % for list of classes
end  % for list of K

% Random tests: ----------------------------------------------------------------
fprintf('\n== Random tests (%g sec):\n', randDelay);

iniTime = cputime;
nTest   = 0;
while cputime - iniTime < randDelay
   K  = round(rand * 8);
   X  = 127 * rand(floor(K + rand * 10 + 1), 1);  % < 127 for INT8
   Y1 = VChooseK(X, K);
   Y2 = nchoosek(X, K);
   if not(isequal(Y1, Y2)) && ~(isempty(Y1) && isempty(Y2))
      error(ErrID, 'Failed for random test data.');
   end
   nTest = nTest + 1;
end
fprintf('  ok: VChooseK equals NCHOOSEK for %d random tests arrays.\n', nTest);

% Invalid input: ---------------------------------------------------------------
fprintf('\n== Check rejection of bad input:\n');
tooLazy = false;

try  % Bad type of empty input:
   dummy   = VChooseK({}, 1);  %#ok<*NASGU>
   tooLazy = true;
catch
   disp(['  ok: {} rejected: ', LF, '      ', strrep(lasterr, LF, '; ')]);
end
if tooLazy
   error(ErrID, '{} not rejected.');
end

try  % Bad type of input:
   dummy   = VChooseK({1, 2}, 1);
   tooLazy = true;
catch
   disp(['  ok: {1, 2} rejected: ', LF, '      ', ...
      strrep(lasterr, LF, '; ')]);
end
if tooLazy
   error(ErrID, '{1, 2} not rejected.');
end

try  % Output too large:
   dummy   = VChooseK(zeros(1, 100), 15);
   tooLazy = true;
catch
   disp(['  ok: Too large input rejected: ', LF, '      ', ...
      strrep(lasterr, LF, '; ')]);
end
if tooLazy
   error(ErrID, 'Too large input not rejected.');
end

try  % Too few inputs:
   dummy   = VChooseK(1);
   tooLazy = true;
catch
   disp(['  ok: No input rejected: ', LF, '      ', ...
      strrep(lasterr, LF, '; ')]);
end
if tooLazy
   error(ErrID, 'No input not rejected.');
end

try  % Too many inputs:
   dummy   = VChooseK(1:2, 3, 4);
   tooLazy = true;
catch
   disp(['  ok: 3 inputs rejected: ', LF, '      ', ...
      strrep(lasterr, LF, '; ')]);
end
if tooLazy
   error(ErrID, '3 inputs not rejected.');
end

% Further ideas for senseful invalid inputs?

% Speed: -----------------------------------------------------------------------
if doSpeed
   fprintf('\n== Speed test (test time: %g sec):\n', SpeedTime);
else
   fprintf('\n== Speed test (test time: %g sec - may be inaccurate):\n', ...
      SpeedTime);
end

% List of function to compare with:
FuncList = { ...
      'NCHOOSEK',    'Matlab toolbox'; ...
      'NCHOOSE2',    'Jos van der Geest, v 2.1 Jun-2008'; ...
      'NChoose2q',   'Jan Simon, Mex, 17-Dec-2009'; ...
      'COMBINATOR',  'Matt Fig, 30-May-2009'; ...
      'VChooseK_M', 'Jan Simon, Matlab version for testing, 21-Dec-2009'; ...
      'VChooseK',   'Jan Simon, Mex, 21-Dec-2009'};

% Loop over data size and input classes:
if doSpeed
   DataClass        = {'int8', 'int16', 'single', 'double'};
   DataClass_sizeof = [1, 2, 4, 8];  % In Bytes
else  % Minimal test only:
   DataClass        = {'double'};
   DataClass_sizeof = 8;  % In Bytes
   
   % Test only NCHOOSEK and VChooseK:
   tmpInd   = ismember(FuncList(:, 1), {'NCHOOSEK', 'VChooseK'});
   FuncList = FuncList(tmpInd, :);
end

% Search for these tools:
Found = true(size(FuncList, 1), 1);
for iFunc = 1:size(FuncList, 1)
   Found(iFunc) = ~isempty(which(FuncList{iFunc, 1}));
end
FuncList = FuncList(Found, :);

% Remove local Matlab version from list of functions if wanted:
if ~testLocalMVersion
   tmpInd              = strcmp(FuncList(:, 1), 'VChooseK_M');
   FuncList(tmpInd, :) = [];
end

% Show function names with authors:
FuncList_T = transpose(FuncList);
fprintf('%-12s:  %s\n', FuncList_T{:});

% For 512 MB, the 8000x1 array needs slow disk caching often. Therefore I
% recommend to include it only, if you have more than 1GB RAM.
DataList = { ...
      2, [5, 10, 100, 1000, 2000, 4000]; ... % , 8000];
      3, [5, 10, 20, 100, 200, 400]; ...
      4, [5, 10, 20, 40, 80, 160]; ...
      5, [5, 10, 20, 40, 80]; ...
      6, [10, 20, 40]};

TooSlow = sprintf(' %6s', 'slow');
Crashed = sprintf(' %6s', 'crash');
Unknown = sprintf(' %6s', '?');

fprintf('\n"slow": omitted, because previous input length took > 1 second\n');
fprintf('For valid tests > 750MB free RAM is needed!\n');

for iData = 1:size(DataList, 1)
   K       = DataList{iData, 1};
   DataLen = DataList{iData, 2};
   
   fprintf('\nK=%d, input length:    ', K);
   fprintf('%7d', DataLen);
   fprintf('\n');
   for iFunc = 1:length(FuncList)
      aFunc = FuncList{iFunc, 1};
      
      % NCHOOSE2 works with K==2 only:
      if K ~= 2 && (strcmpi(aFunc, 'NCHOOSE2') || strcmpi(aFunc, 'NChoose2q'))
         continue;
      end
      
      for iClass = 1:length(DataClass)
         fprintf('  %-12s %-6s ', aFunc, DataClass{iClass});
         NPerSec = Inf;
         
         for aDataLen = DataLen
            if NPerSec < 1  % Give up, if the former loop was too slow:
               fprintf(TooSlow);
               continue;
            end
            
            % Create vector for specific class:
            aData = my_cast(mod(1:aDataLen, 127), DataClass{iClass});
            N     = 0;
            try
               switch aFunc
                  case 'NCHOOSEK'
                     % Matlab toolbox
                     iTime = cputime;
                     while cputime - iTime < SpeedTime
                        S = nchoosek(aData, K);  clear('S');
                        N = N + 1;
                     end
                     NPerSec = N / (cputime - iTime);  % Loops per second
                     PrintLoop(NPerSec);
                     
                  case 'NCHOOSE2'
                     % FEX: Jos van der Geest, version 2.1 (jun 2008)
                     iTime = cputime;
                     while cputime - iTime < SpeedTime
                        S = nchoose2(aData);  clear('S');
                        N = N + 1;
                     end
                     NPerSec = N / (cputime - iTime);  % Loops per second
                     PrintLoop(NPerSec);
                     
                  case 'NChoose2q'
                     % Mex tuned for 2 elements (was temporarily on the FEX):
                     iTime = cputime;
                     while cputime - iTime < SpeedTime
                        S = NChoose2q(aData);  clear('S');
                        N = N + 1;
                     end
                     NPerSec = N / (cputime - iTime);  % Loops per second
                     PrintLoop(NPerSec);
                     
                  case 'COMBINATOR'
                     % Matt Fig, FEX, 5/30/2009
                     iTime = cputime;
                     while cputime - iTime < SpeedTime
                        S = combinator(aDataLen, K, 'c');  clear('S');
                        N = N + 1;
                     end
                     NPerSec = N / (cputime - iTime);  % Loops per second
                     PrintLoop(NPerSec);
                     
                  case 'VChooseK'
                     % FEX: Jan Simon (this function is tested!)
                     iTime = cputime;
                     while cputime - iTime < SpeedTime
                        S = VChooseK(aData, K);  clear('S');
                        N = N + 1;
                     end
                     NPerSec = N / (cputime - iTime);  % Loops per second
                     PrintLoop(NPerSec);
                     
                  case 'VChooseK_M'
                     % Matlab implementation as proof of concept:
                     iTime = cputime;
                     while cputime - iTime < SpeedTime
                        S = VChooseK_M(aData, K);  clear('S');
                        N = N + 1;
                     end
                     NPerSec = N / (cputime - iTime);  % Loops per second
                     PrintLoop(NPerSec);
                     
                  otherwise
                     error(ErrID, 'Unknown function.');
               end
               
            catch
               % Most likely the called function is not compatible to the used
               % Matlab version (happened for 6.5 only):
               fprintf(Crashed);
            end
            drawnow;  % Give CTRL-C a chance!
         end  % for aDataLen
         fprintf('  loops/sec\n');
         
         % Extra pause, because the hard disk could be used as virtual memory:
         if NPerSec < 0.1
            pause(5);
         end
      end  % for DataLen
   end  % for iFunc
end  % for iData: DataList

fprintf('\n==> VChooseK seems to work fine.\n');

return;

% ******************************************************************************
function PrintLoop(N)
if N > 10
   fprintf(' %6.0f', N);
else
   % fprintf(' %6.1f', N);
   fprintf(' %6.2f', N);
end

return;

% ******************************************************************************
function A = my_cast(A, ClassName)
% Simulate CAST for Matlab 6
A = feval(ClassName, A);

return;

% ******************************************************************************
function Y = VChooseK_M(X, K)
% VChooseK - Combinations of K elements [MEX]
% See VChooseK for help.
% Although this was thought to debug the C-Mex algorithm only, it is e.g. 60
% times faster then NCHOOSEK(1000, 5)!
% In Matlab 7.8 the function run fastest for DOUBLEs.
%
% Tested: Matlab 6.5, 7.7, 7.8, WinXP, [UnitTest]
%         Compilers: BCC5.5, LCC2.4/3.8, Open Watcom 1.8
% Author: Jan Simon, Heidelberg, (C) 2009-2010 J@n-Simon.De
% License: BSD (use/copy/modify on own risk, but mention author)

% ==============================================================================
% Slower M-version as proof of concept - no input checks here.
X  = X(:);         % Give X a column shape
nX = numel(X);
if nX <= K         % Check pathological input
   if nX == K
      Y = transpose(X);
   else
      Y = [];
   end
   return;
end

switch K
   case 1
      Y = X;
      
   case 2   % 100 times faster than NCHOOSEK(1:1000, 2)!
      Y  = zeros(nX * (nX - 1) / 2, 2);
      a  = 1;
      for k = 1:nX - 1
         b         = a + nX - k - 1;
         Y(a:b, 1) = X(k);
         Y(a:b, 2) = X(k + 1:nX);
         a         = b + 1;
      end
      
   case 3  % Just insert another loop:
      Y  = zeros(nX * (nX - 1) * (nX - 2) / 6, 3);
      a  = 1;
      for k1 = 1:nX - 2
         for k2 = k1 + 1:nX - 1
            b         = a + nX - k2 - 1;
            Y(a:b, 1) = X(k1);
            Y(a:b, 2) = X(k2);
            Y(a:b, 3) = X(k2 + 1:nX);
            a         = b + 1;
         end
      end
      
   case 4   % Insert another loop again:
      % 4.4 times faster than Matlab's NCHOOSEK, but the MEX is 240 times faster
      % for (1:20).
      nY = round(((((nX * (nX - 1) / 2) * (nX - 2)) / 3) * (nX - 3)) / 4);
      Y  = zeros(nY, 4);
      a  = 1;
      for k1 = 1:nX - 3
         for k2 = k1 + 1:nX - 2
            for k3 = k2 + 1:nX - 1
               b         = a + nX - k3 - 1;
               Y(a:b, 1) = X(k1);
               Y(a:b, 2) = X(k2);
               Y(a:b, 3) = X(k3);
               Y(a:b, 4) = X(k3 + 1:nX);
               a         = b + 1;
            end
         end
      end
      
   otherwise
      % It is trivial to expand this more and more by inserting loops over
      % k4, k5, k6... The a general algorithm is more challenging:
      
      % Use a loop to calculate the length ofthe output. All intermediate
      % values are integers to reduce rounding problems.
      d  = nX - K;
      nY = d + 1;
      for i = 2:K
         nY = nY + (nY * d) / i;
      end
      
      % Get memory for output:
      Y = zeros(round(nY), K);
      
      % The algorithm is designed in C-style, because it was used for debug the
      % C-Mex function. Just the loop over the K.th column is vectorized. With
      % using the power of Matlab, a much faster algorithm would be possible
      % (see e.g. COMBINATOR of Matt Fig, FEX: 24325). But this is just a dummy
      % for the much faster MEX!
      Index = 1:K;                       % Current index
      Limit = [(nX - K + 1):nX - 1, 0];  % Last index not checked here!
      a     = 1;
      while 1
         b = a + nX - Index(K);          % Write index for last column
         for i = 1:(K - 1)               % Write the left K-1 columns
            Y(a:b, i) = X(Index(i));
         end
         Y(a:b, K) = X(Index(K):nX);     % Write the K.th column
         a         = b + 1;              % Move the write pointer
         
         % Search the last column index, which is not exhausted:
         newLoop = find(Index < Limit);
         if isempty(newLoop)             % All columns are filled:
            return;                      % Ready!
         end
         newLoop = newLoop(length(newLoop));
         
         % Fill new Index with new value encreasing by 1:
         Index(newLoop:K) = Index(newLoop) + (1:(K - newLoop + 1));
      end
end  % switch K

return;
