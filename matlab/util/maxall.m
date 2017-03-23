function allmax = maxall(I)
% MAXALL Maximum value in N-dimensional array.
%   MAXALL(I) returns the largest value in the array I.
%
%   For example:
%   MaxGreenVal = maxall(Irgb(:,:,2));
%
%   Another example:
%   MaxDims2and5 = maxall(data(:,[2 5]));
%
%   See also STRAIGHTEN, MINALL, MEANALL, ALLL
%
%JC

allmax = max(I(:));
