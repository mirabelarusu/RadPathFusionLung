function allsum = sumall(I)
% SUMALL Sum of values in N-dimensional array.
%   SUMALL(I) returns the sum of values in the array I.
%
%   For example:
%   SumGreenVal = sumall(Irgb(:,:,2));
%
%   Another example:
%   SumDims2and5 = sumall(data(:,[2 5]));
%
%   See also STRAIGHTEN, MINALL, MEANALL, ALLL
%
%JC

allsum = sum(I(:));
