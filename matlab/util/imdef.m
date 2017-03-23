function [Ifwd,M,invM,Irev,cpointsIout] = imdef(params,I,method,fillwith,preM,coordshift,cpointsI,framesize)
% IMDEF Deform an image or image set with affine transformation, with 
%   origin at center of image or user specified location, returning an
%   image of equal size as the input.
%
%   IMDEF(PARAMS,I) deforms image I subject to deformations specified
%   by the transformation parameters PARAMS. PARAMS is a 1D array with the
%   elements corresponding to:
%           [rotation, trans_x, trans_y, scale_x, scale_y]
%   Either 1, 3, or 5 parameters may be specified.
%
%   IMDEF(M,I) performs the transformation specified by the linear 
%   coordinate transformation matrix M. M is 3-by-3 or 4-by-4 homogeneous
%   and representing an in-plane (x-y) transformation.
%
%   Optional INPUT arguments (3 through 6):
%
%   IFWD = IMDEF(...,METHOD) uses the specified interpolation METHOD.
%   Possible values are 'nearest', 'linear', and 'cubic'. Both 'linear' and
%   'cubic' can return values not present in I, and 'cubic' can return
%   values outside the range of values in I. Default METHOD is 'nearest'. 
%   Specify [] for default.
%
%   IFWD = IMDEF(...,METHOD,FILLWITH) uses background value specified by
%   scalar FILLWITH. Default is zero (0). Specify [] for default.
%
%   IFWD = IMDEF(...,METHOD,FILLWITH,PREM) first applies the 4-by-4
%   homogeneous transformation matrix PREM. Default is eye(4). Specify [] 
%   for default.
%
%   IFWD = IMDEF(...,METHOD,FILLWITH,PREM,COORDSHIFT) uses COORDSHIFT to 
%   specify an offset from the frame center about which deformations are 
%   performed. Positive values impose a shift toward the top and left.
%   Default COORDSHIFT is [0; 0]. !![Y; X]!! Specify [] for default.
%
%   TIP: To rotate about the image origin (top left corner), use a value of
%   COORDSHIFT = [size(I,1)/2 size(I,2)/2] to affect a shift from the image
%   center to the top left corner.
%
%   Optional OUPUT arguments:
%
%   [IFWD, M, INVM] = IMDEF(...) returns the transformation matrix M and 
%   the inverse INVM.
%
%   [IFWD, M, INVM, IREV] = IMDEF(...) the result of IREV applied to I.
%
%   [IFWD, M, INVM, IREV, CPOUT] = IMDEF(...,CPIN) transforms the control 
%   point coordinates in the 7th input argument CPIN, returning the new 
%   coorinates in CPOUT. CPIN is a 2-by-N array (i.e. [x; y]).
%
%   [IFWD, M, INVM, IREV, CPOUT] = IMDEF(...,CPIN,FRAMESIZE) transforms 
%   CPIN in an image frame of size FRAMESIZE, if I is empty. If I is not
%   empty, FRAMESIZE is *ignored*, a warning is issued, and CPIN is 
%   transformed about the center of the image frame.
%
%   NOTES:
%   -Fixed output image frame size. The output image is the same size
%    as the input image. For an equal or larger output, see IMDEFFULL.
%   -If I is not a float, it is converted to double.
%
%   Example of nearest neighbor interpolation:
%   I=imread('cameraman.tif'); I=double(I);
%   cornerbraces=sub2ind(size(I),...
%   [1 1 2 1 1 2 size(I,1) size(I,1) size(I,1)-1 size(I,1) size(I,1) size(I,1)-1],...
%   [1 2 1 size(I,2) size(I,2)-1 size(I,2) 1 2 1 size(I,2) size(I,2)-1 size(I,2)]);
%   I(cornerbraces)=max(I(:));
%   [Ifwd,M,invM,Irev]=imdef([.1 10 -6 1.1 .9],I,'nearest');
%   figure; imagesc(I); colormap gray; axis off equal
%   figure; imagesc(Ifwd); colormap gray; axis off equal
%   figure; imagesc(Irev); colormap gray; axis off equal
%   isequal(size(I),size(Ifwd))   % true
%
%   Example of linear interpolation (load image as above):
%   [Ifwd,M,invM,Irev]=imdef([.1 10 -6 1.1 .9],I,'linear');
%   figure; imagesc(I); colormap gray; axis off equal
%   figure; imagesc(Ifwd); colormap gray; axis off equal
%   figure; imagesc(Irev); colormap gray; axis off equal
%   isequal(size(I),size(Ifwd))   % true
%
%   Class support for I:  
%      float: double, single
%   
%   See also ROTMAT, IMDEFFULL, INTERP2.
%
%JC

% This function uses s 4-by-4 homogeneous transformation matrix to
% implement a 2D deformation instead of using a 3-by-3 matrix. This is done
% in order to facilitate code sharing with a 3D transformation framework.

% Ensure I is a float
if ~isa(I,'float'),
    I=double(I);
    warning('imdef:ConvertingToDouble','Converting I to double datatype.');
end
Iclass=class(I);

% Check image size
[nrows, ncols, planes]=size(I); % note that planes may be >size(I,3). fine.
% [nrows,ncols,spatialplanes,attributeplanes]=size(I);
% planes=spatialplanes*attributeplanes;

% Parse Inputs for: coordshift, preM, fillwith, and method
% correcting=true;
if nargin<6 || isempty(coordshift), % no coordshift
    coordshift=zeros(2,1); % zeros(2,1,Iclass); % y, x
    % correcting=false;
    if nargin<5 || isempty(preM), % no preM
        preM=eye(4); % eye(4,Iclass);
        if nargin<4 || isempty(fillwith), %~exist('fillwith','var'),
            fillwith=zeros(1,1,planes,Iclass); %cast(0,class(I));
            if nargin<3 || isempty(method),
                method='nearest';
            end
        end
    end
end

if isempty(coordshift), coordshift=zeros(2,1); %correcting=false;
    end
if isempty(fillwith), fillwith=zeros(1,1,planes,Iclass); end
if isempty(preM), preM=eye(4); end
if isempty(method), method='nearest'; end

% if isa(I,'single'),
%     coordshift=single(coordshift);
%     preM=single(preM);
%     fillwith=single(fillwith);
%     params=single(params);
%     nrows=single(nrows);
%     ncols=single(ncols);
% end
npts=nrows*ncols;

% Check transformation parameter input
if min(size(params,1),size(params,2))==1, % got parameters for rotation, etc. instead of matrix
    if length(params)<5,
        sx=1;sy=1;
        if length(params)<3,
            dx=0;dy=0;
        else
            dx=params(2);
            dy=params(3);
        end
    elseif length(params)>6,
        error('Too many parameters.  Specify 1, 3, 5 or 6 only.');
    else
        sx=params(4);
        sy=params(5);
        dx=params(2);
        dy=params(3);
    end
    theta=params(1);
    if numel(params)==6, theta2=params(6);
    else theta2=0; end

    % Homogeneous translation matrix
    transM=eye(4);
    transM(:,4)=[dx dy 0 1]';

    % Homogeneous scaling matrix
    scaleM=diag([sx sy 1 1]);

    % Homogeneous rotation matrix, about z-axis
    rotM=[cos(theta) sin(theta) 0 0;
        -sin(theta) cos(theta) 0 0;
        0          0           1 0;
        0          0           0 1];
    rotM2=[cos(theta2) sin(theta2) 0 0;
        -sin(theta2) cos(theta2) 0 0;
        0          0           1 0;
        0          0           0 1];

    % Rotate, Scale, (Rotate,) then Translate
    M=transM*rotM2*scaleM*rotM*preM;

elseif size(params,1)==3,  % got 3-by-3 transformation matrix, make homogeneous
    M=eye(4);
    M(1:2,1:2)=params(1:2,1:2);
    M(1:2,4)=params(1:2,3);
    M=M*preM;
elseif size(params,1)==4,  % got 4-by-4 homogeneous transformation matrix
    M=params;
    M=M*preM;
else
    error('Matrix not 3x3 or 4x4 homogneous.');
end

% Inverse operation so every target pixel gets a value, speedy interp2 assignment
% M=cast(M,Iclass); % should not need this
invM=inv(M);

% check if image was [] (i.e. only control points input)
if npts==0,
    nrows=framesize(1);
    ncols=framesize(2);
elseif exist('framesize','var') && ~isempty(framesize),
    warning('imdef:framesizeIgnored','Image not empty, framesize variable IGNORED!');
end

% Center of image, adjusted center of image
centpt=[(ncols+1)/2-coordshift(2); (nrows+1)/2-coordshift(1)];  % x, y
% if correcting,
%     impcentpt=[(ncols/2); (nrows/2)];
%     fprintf('Implict center at (X,Y) = (%g,%g).\n',impcentpt(1),impcentpt(2));
%     fprintf('Transforming about center at (X,Y) = (%g,%g).\n',centpt(1),centpt(2));
% end

% Want control points?
cpointsIout=[];
if nargin>6 && nargout>4,
    % Make CPs homogeneous
    if size(cpointsI,1)==2,
        cpointsI=[cpointsI; ones(2,size(cpointsI,2))];
    elseif size(cpointsI,1)==3,
        cpointsI=[cpointsI; ones(1,size(cpointsI,2))];
    end
    
    % Move origin to image center for sensible rotation
    cpointsI(1,:)=cpointsI(1,:)-centpt(1); % x
    cpointsI(2,:)=cpointsI(2,:)-centpt(2); % y
    % move control points
    cpointsIout=M*cpointsI;
    % move grid back to matlab coordinates
    cpointsIout(1,:)=cpointsIout(1,:)+centpt(1); % x
    cpointsIout(2,:)=cpointsIout(2,:)+centpt(2); % y
    % remove homogeneous
    cpointsIout(4,:)=[];
end

% Don't work on image?
if isempty(I),
    Ifwd=[]; Irev=[];
    return;
end

% Produce coordinate grid
[xpts,ypts]=meshgrid(1:ncols,1:nrows);

% Make homogeneous coordinates vector
% homog_points=[xpts(:)'; ypts(:)'; ones(2,npts)];
homog_points_centered=ones(4,npts); % ones(4,npts,Iclass);

% Put image center at (0,0) for sensible rotation
homog_points_centered(1,:)=xpts(:)-centpt(1);
homog_points_centered(2,:)=ypts(:)-centpt(2);
clear xpts ypts

% Apply the transform
% moved_points=invM*homog_points_centered;
moved_points=M\homog_points_centered;
if nargout<4, clear homog_points_centered; end

% Move the grid back to matlab coordinates
newx=moved_points(1,:)+centpt(1);
newy=moved_points(2,:)+centpt(2);
clear moved_points

% Do forward transformation of image
Ifwd=zeros(size(I),Iclass);
if numel(fillwith)~=planes, fillwith=repmat(fillwith,planes,1); end
if numel(fillwith)~=planes, error('fillwith input wrong size'); end
for iz=1:planes,
%     Inew(:,:,iz)=griddata(reshape(newx,[nrows ncols]),reshape(newy,[nrows ncols]),I(:,:,iz),xpts,ypts,'nearest');
    if strfind(method,'linear'),
        Ifwd(:,:,iz)=linterp2mex(I(:,:,iz),reshape(newx,[nrows ncols]), ... 
        reshape(newy,[nrows ncols]),fillwith(iz));
    else
        Ifwd(:,:,iz)=interp2(I(:,:,iz),reshape(newx,[nrows ncols]), ... 
        reshape(newy,[nrows ncols]),method,fillwith(iz));
    end
end
% Ifwd(isnan(Ifwd))=fillwith; % should not be necessary with EXTRAPVAL=fillwith
for iz=1:planes,
    wherenans=any(isnan(Ifwd),3);
    Ifwd(wherenans)=fillwith(iz); % should not be necessary with EXTRAPVAL=fillwith
end

% Do reverse transformation of image
if nargout>3,
    forward_points=M*homog_points_centered;
    forward_points(1,:)=forward_points(1,:)+centpt(1);
    forward_points(2,:)=forward_points(2,:)+centpt(2);
    
    Irev=zeros(size(I),Iclass);
    for iz=1:planes,
        Irev(:,:,iz)=interp2(I(:,:,iz),reshape(forward_points(1,:),[nrows ncols]), ... 
            reshape(forward_points(2,:),[nrows ncols]),method,fillwith(iz));
    end
    % Irev(isnan(Irev))=fillwith; % should not be necessary with EXTRAPVAL=fillwith
    for iz=1:planes,
        wherenans=any(isnan(Irev),3);
        Irev(wherenans)=fillwith(iz); % should not be necessary with EXTRAPVAL=fillwith
    end
end
