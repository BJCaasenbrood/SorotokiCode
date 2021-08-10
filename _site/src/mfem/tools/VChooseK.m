function Y = VChooseK(X, K)
% VChooseK - Combinations of K elements [MEX]
% VChooseK(V, K) creates a matrix, which rows are all combinations of 
% choosing K elements of the vector V without order and without repititions.
%
% INPUT:
%   V: Array of class DOUBLE, SINGLE, (U)INT8/16/32/64, LOGICAL, CHAR.
%      The shape of the array does not matter.
%   K: Number of elements to choose. Scalar double with integer number.
%
% OUTPUT:
%   Y: Matrix of size [N!/K!(N-K)!, K] with N is the number of elements of V.
%      Y has the same class as the input V. If V has less than K elements, a
%      warning appears and the empty matrix is replied. As usual in Matlab,
%      empty inputs yield to empty Y.
%      A warning appears, if the output would exceed 500MB, and an error for
%      1GB. Both limits can be adjusted according to the available RAM in the
%      C-Mex source.
%
% NOTE: The output equals Matlab's NCHOOSEK, except that NCHOOSEK replies the
% number of combinations for a scalar V:
%   VChooseK(-1, 1) replies [-1]: One element taken out of a set of length one.
%   NCHOOSEK(-1, 1) fails at calculating N!/K!(N-K)!.
%
% EXAMPLES:
%   Choose 2 elements from [1,2,3,4]:
%     VChooseK(1:4, 2)
%     ==> [1,2; 1,3; 1,4; 2,3; 2,4; 3,4]
%   For speed cast the input to integer types if possible:
%     Y = double(VChooseK(int16(1:1000), 2);
%   is faster than:
%     Y = VChooseK(1:1000, 2);
%   To get the combinations of cell arrays, use the combinations of the index:
%     C  = {'a', 'b', 'c', 'd'};
%     C2 = C(VChooseK(1:4, 2))
%     ==> C2 = {'a', 'b'; 'a', 'c'; 'a', 'd'; 'b', 'c'; 'b', 'd'; 'c', 'd'}
%
% COMPILE:
%   mex -O VChooseK.c
% Compatibility to 64-bit machines is assumed, but cannot be tested currently.
% On Linux the C99 comments must be considered (thanks Sebastiaan Breedveld):
%   mex CFLAGS="\$CFLAGS -std=C99" -O VChooseK.c
% Please run the unit-test TestVChooseK after compiling!
%
% INSPIRATION:
% Jos van der Geest has published an efficient NCHOOSE2, which is much faster
% than Matlab's NCHOOSEK:
%   http://www.mathworks.com/matlabcentral/fileexchange/20144
% I was curious, if a MEX version is much faster. At first I've created
% NChoose2q, which is 5 times faster than Jos' Matlab implementation. But the
% algorithm was so simple in C and did not consume temporary memory, that I've
% expanded it to K=3,4,5 soon. The algorithm for free K was more challenging,
% but it needs just 3*K pointers for temporary memory. Therefore it is much
% faster than Matlab's NCHOOSEK (Matlab 6.5, 1.5GHz Pentium-M, 512 MB):
%                    NCHOOSEK:   VChooseK:
%   (1:50, 5):       45.30 sec   0.32 sec      => 141 times faster
%   (1:20, 16):       0.52 sec   0.0024 sec    => 216 times faster
%   (int8(1:40), 4):  0.64 sec   0.001562 sec  => 410 times faster
%
% There is a Matlab implementation in TestVChooseK, which might be helpful to
% understand the applied algorithm.
%
% Related FEX related publications:
% (http://www.mathworks.com/matlabcentral/fileexchange/<number>)
%   NCTWO (Simone Scaringi): 20110
%   COMBINATOR (Matt Fig) [I prefer this for general tasks!]: 24325
%   PICK (Stefan Stoll) calls NCHOOSEK for unordered combinations: 12724
%   COMBINATIONS (Gautam Vallabha): 23080
%
% Tested: Matlab 6.5, 7.7, 7.8, WinXP, [UnitTest]
%         Compilers: BCC5.5, LCC2.4/3.8, Open Watcom 1.8
% Author: Jan Simon, Heidelberg, (C) 2009 matlab.THIS_YEAR(a)nMINUSsimonDOTde
% License: BSD (use/copy/modify on own risk, but mention author)

% $JRev: R0g V:005 Sum:yQLycAwkArU7 Date:24-Dec-2008 02:59:38 $
% $File: User\JSim\Published\VChooseK\VChooseK.m $
% History:
% 001: 17-Dec-2009 08:10, NChoose2q: MEX function for NCHOOSEK(X, 2).
% 004: 21-Dec-2009 13:49, VChooseK: Expanded for free K.

% There is a Matlab implementation in TestVChooseK.
