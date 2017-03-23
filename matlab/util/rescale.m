function [Iout,drange,dmin] = rescale(I)
% RESCALE Rescale data into the range [0,1].
%   RESCALE(I) rescales the array I so that all elements fall in the range
%   [0,1]. The output is double precision.
%
%   See also RESCALE_RANGE.
%
%JC

% Convert input to double precision
% if ~isa(I,'double'),
%     I=double(I);
% end
if ~isfloat(I),
    I=double(I);
elseif isa(I,'single'),
    warning('rescale:gotSingleData','Input data is single precision.  Results may be inaccurate.');
end

% Make sure the data can be rescaled with the current machine precision
dmin=min(I(:));
drange=max(I(:))-dmin;
if drange > eps,
    % Iout = (I- min(I)) / range(I)
    Iout = (I-dmin)/drange;
else
    Iout=I;
end
