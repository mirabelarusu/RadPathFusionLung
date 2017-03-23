function Ibound = impad(I,bsize,padwith)
% IMPAD Pad an image.
%   IMPAD(I,BSIZE,PADWITH) returns the image I padded to size BSIZE. The
%   scalar value specified by PADWITH is used as background.
%
%   For example:
%   xg = imread('cameraman.tif');
%   Ibound = impad(xg,[320 320],cast(0,class(xg)));
%   imshow(Ibound)
%
%JC

brows=bsize(1);
bcols=bsize(2);
Isize=size(I);
Idims=ndims(I);
nrows=Isize(1); ncols=Isize(2);
if Idims>2, otherdims=Isize(3:Idims); else otherdims=1; end
nhighdims=prod(otherdims);

% Default to pad with zero of same class as I
if nargin<3,
    padwith=cast(0,class(I));
end

ydiff=brows-nrows;
xdiff=bcols-ncols;

if ydiff<0 || xdiff<0, error('Pad size smaller than original size.'); end

% Imiddlex=ncols/2;
% Imiddley=nrows/2;

yborder=floor(ydiff/2);
xborder=floor(xdiff/2);

% Ibound=padwith*ones(bsize(1),bsize(2),size(I,3));
if numel(padwith)==size(I,3),
    padwith=permute(padwith(:),[3 2 1]);
    Ibound=repmat(padwith,[bsize(1) bsize(2) 1]); % possibly >3D
elseif numel(padwith)==1,
    Ibound=repmat(padwith,[bsize(1) bsize(2) otherdims]); % possibly >3D
else
    error('Dimensions of padwith do not agree with image depth.');
end

for i=1:nhighdims,
    Ibound((1:nrows)+yborder,(1:ncols)+xborder,i) = I(:,:,i);
end

% for i=1:size(I,3),
%     Ibound((1:brows-ydiff)+yborder,(1:bcols-xdiff)+xborder,i) = ... 
%         I((1:nrows-ydiff)+yborder,(1:ncols-xdiff)+xborder,i);
% end
