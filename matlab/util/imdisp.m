function [varargout]=imdisp(varargin)
% IMDISP Display and image with reasonable coloring and equal aspect ratio.
%   IMDISP(I) displays image I.
%
%   IMDISP(I,CMAP) displays indexed image I with colormap CMAP.
%
%   IMDISP(HA,I) displays image I on axes with handle HA.
%
%   IMDISP(HA,I,CMAP) displays index image I with colormap CMAP on axes HA.
%
%   This is a good replacement for imshow if you don't have it.
%
%   Note: also works for RGB images, but one of two things must be true:
%   (1) double precision must be scaled on [0,1]
%   (2) integer must be scaled on [0, 255]
%
%JC

% drawnow; % in case you just made a figure and/or axes
gotaxhand=false;
if numel(varargin{1})==1 && ishandle(varargin{1}),
    gotaxhand=true;
elseif numel(varargin{1})==1,
    drawnow;
    if ishandle(varargin{1}), gotaxhand=true; end
end

if nargin>1 && gotaxhand,
    % ax=varargin{1};
    % I=varargin{2};
    if isempty(varargin{2}), error('Empty image.'); end
    if nargin>2,
        hi=image(varargin{2},'CDataMapping','direct','Parent',varargin{1});
        set(gcf, 'Colormap', varargin{3});
        if size(varargin{3},1)>256 && ispc, set(gcf,'Renderer','Zbuffer'); end
    else
        hi=image(varargin{2},'CDataMapping','scaled','Parent',varargin{1});
        colormap gray;
    end
elseif nargin==2,
    hi=image(varargin{1},'CDataMapping','direct');
    set(gcf, 'Colormap', varargin{2});
    if size(varargin{2},1)>256 && ispc, set(gcf,'Renderer','Zbuffer'); end
elseif nargin==1,
    % I=varargin{1};
    if isempty(varargin{1}), error('Empty image.'); end
    hi=imagesc(varargin{1});
    colormap gray; 
else
    error('Wrong inputs');
end
if nargout==1,
    varargout{1}=hi;
end

axis off equal tight
