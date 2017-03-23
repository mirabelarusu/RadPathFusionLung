function Icrop = imshave(I,bsize)
% IMSHAVE Crop an image evenly on all sizes.
%   IMSHAVE(I,BSIZE) returns the image I cropped to size BSIZE.
%
%   For example:
%   xg = imread('cameraman.tif');
%   Icrop = imshave(xg,[200 200]);
%   imshow(Icrop)
%
%JC

brows=bsize(1);
bcols=bsize(2);
nrows=size(I,1);
ncols=size(I,2);
if any([brows bcols]>[nrows ncols]),
    error('New size [%g,%g] must be smaller than original image size [%g,%g].',brows,bcols,nrows,ncols);
end

if ndims(I)>8,
    warning('imshave:HighDimsCollapsed','Dimensions above 8 will be collapsed into dimension 8.');
end

ydiff=nrows-brows;
xdiff=ncols-bcols;

yborder=floor(ydiff/2);
yextra=ceil(rem(ydiff/2,1));
xborder=floor(xdiff/2);
xextra=ceil(rem(xdiff/2,1));

% Icrop=zeros([bsize(1),bsize(2),size(I,3)]);
Icrop=I(yborder+1:end-yborder-yextra,xborder+1:end-xborder-xextra,:,:,:,:,:,:);
