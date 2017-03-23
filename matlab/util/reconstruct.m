function [reconstruction, reconstructionR, reconstructionG, reconstructionB, mask_4D_reg, rigid_params] ...
    = reconstruct ( path_hist, path_out, varargin )
% reconstruct - 3D reconstruction from 2D stacks of images
%
%   path_hist - where the input histology images are located on the HD
%   path_out  - where to output register images should be wrote out
%   path_masks - where masks are located
%       note: masks have to have the exact same name as the histology files
%   prefix_masks - all important marks are located: can be folders
%   
%   id_roi which mask should be used to drive the registration (used to crop the img)
%
%   newSize - size that fits all images
%   magni - how much to downsample (=5 for 5 time downsampling)
%   order - data is read in alfabetic order, but that might not be what one desires
%       default should be order = 1:4;
%       example: order = [3,1,2,4];
%   voxel_size - Final reconstruction voxel size: default 1,1,1 should reflect all
%       scaling, magnification, shrinking etc ...
%   flipZ - should the data be flipped around Z?
%   PaddingOnZ - how much padding should be added when experting the mha
%   weight_entire_sample - weight of the entire sample
%   ignore_padding_target: How many slices to ignore in target (number of 
%                  slices the histology volume was padded with)
%   score_forward: should the forward slices be considered in the scoring
%
%   function
%   id_sample
ip=inputParser;
ip.addParamValue('path_masks', '', @(x)ischar(x));
ip.addParamValue('prefix_masks', '', @(x)isvector(x));
ip.addParamValue('path_target', '', @(x)ischar(x));
%parameters
ip.addParamValue('hist_prefix', 'image_', @(x)ischar(x));
ip.addParamValue('mask_prefix', 'mask_', @(x)ischar(x));
ip.addParamValue('id_roi', 0, @(x)isnumeric(x) && numel(x)==1);
ip.addParamValue('id_lobe', 1, @(x)isnumeric(x) && numel(x)==1);
ip.addParamValue('new_size', [0,0], @(x)isnumeric(x) && numel(x)==2);
ip.addParamValue('magni', 1, @(x)isnumeric(x) && numel(x)==1);
ip.addParamValue('image_order', [], @(x)isnumeric(x));
ip.addParamValue('weight_entire_sample', 0, @(x)isnumeric(x) && numel(x)==1);
ip.addParamValue('id_sample', 0, @(x)isnumeric(x) && numel(x)==1);
ip.addParamValue('ignored_padding_target', 0, @(x)isnumeric(x) && numel(x)==1);
ip.addParamValue('score_forward', false, @(x)islogical(x) && numel(x)==1);
ip.addParamValue('lobe_ids', {}, @(x)iscell(x));
ip.addParamValue('lobe_scale', 1 , @(x)isnumeric(x) && numel(x)==1);
ip.addParamValue('params', [], @(x)isnumeric(x) && numel(x)>2);
ip.addParamValue('hist_extension', 'tif', @(x)ischar(x));
ip.addParamValue('flip_LR_vol', 0 , @(x)isnumeric(x) && numel(x)==1);
ip.parse(varargin{:});


path_masks = ip.Results.path_masks;
prefix_masks = ip.Results.prefix_masks;
path_target = ip.Results.path_target;
image_order = ip.Results.image_order;
id_roi = ip.Results.id_roi;
id_sample = ip.Results.id_sample;
new_size = ip.Results.new_size;
magni = ip.Results.magni;
weight_entire_sample = ip.Results.weight_entire_sample;
ignored_padding_target = ip.Results.ignored_padding_target;
score_forward  = ip.Results.score_forward;
params_2_refine = ip.Results.params;
lobe_ids = ip.Results.lobe_ids;
lobe_scale = ip.Results.lobe_scale;
hist_prefix = ip.Results.hist_prefix;
mask_prefix = ip.Results.mask_prefix;
hist_extension = ip.Results.hist_extension;
flip_LR_vol = ip.Results.flip_LR_vol;

%hard coded parameters:
N = 64; % number of bins in MI: FIXME: is still hardcoded
K = 4;  % number of considered consecutive slices

pad = 1;
rota = 0;

%%%
%%% Read the input
%%%
[hist_3D, hist_3D_R, hist_3D_G, hist_3D_B, mask_4D, target] = ...
    readReconstructInput (path_hist, hist_prefix, mask_prefix, image_order, new_size, ...
    path_masks, prefix_masks, path_target, N, (weight_entire_sample>0), ...
    ignored_padding_target, hist_extension, flip_LR_vol);

%%%
%%% preprocess the input: rotation and scaling
%%%
hist_3D = preprocessRecontructionInput(hist_3D, magni, rota, 0);
hist_3D_R = preprocessRecontructionInput(hist_3D_R, magni, rota, 0);
hist_3D_G = preprocessRecontructionInput(hist_3D_G, magni, rota, 0);
hist_3D_B = preprocessRecontructionInput(hist_3D_B, magni, rota, 0);
mask_4D = preprocessRecontructionInput(mask_4D, magni, rota, 0);


%%%
%%% execute the reconstruction
%%%

%%% figure out that is the template: mask 
if (id_sample>0 && id_sample<=size(mask_4D,4))
    disp(['Considering both region of interest ', num2str(id_roi), ...
        ' and subregion of interest ', num2str(id_sample), ' with weight ',...
          num2str(weight_entire_sample)]);
    mask = double((mask_4D(:,:,:,id_roi)>0).*(1-weight_entire_sample)+...
        weight_entire_sample*(mask_4D(:,:,:,id_sample)>0));
else
    if (id_roi>0 && id_roi<=size(mask_4D,4))
        disp(['Considering the region of interest ', num2str(id_roi)])
        mask = (mask_4D(:,:,:,id_roi)>0);
    else
        mask = sum(mask_4D,4)>0;
    end
end

if (size(lobe_ids,2)==0)
    lobe_mask(:,:,:,1) = mask;
end

%%% deal with lobes and units
se = strel('disk',1);
if (size(lobe_ids,2)>0) % lobes were provided
    for (i = 1:size(lobe_ids,2))
        ids = lobe_ids{i};
        lobe_mask(:,:,:,i) = imopen(uint8(...
            rescale(mask_4D(:,:,:,id_roi))*lobe_scale)==ids(1),se);
        for (j = 2:size(ids,2))
            lobe_mask(:,:,:,i) = lobe_mask(:,:,:,i)+...
                imopen((uint8(rescale(mask_4D(:,:,:,id_roi))*lobe_scale)==ids(j)),se);
        end
        lobe_mask(:,:,:,i) = lobe_mask(:,:,:,i)>0;
    end
end

for (i = 1:size(lobe_mask,4))
    %%% Create background: considered in the registration, but not
    %%% transformed
    lobe_bg = zeros(size(hist_3D));
    for (j = 1:size(lobe_ids,2))
        if (j<i)
            lobe_bg = uint8(lobe_bg) + uint8(applyTransf( double(hist_3D), ...
                double(lobe_mask(:,:,:,j)), rigid_params(:,:,j), 'linear', 0));
        end
        if (j>i)
            lobe_bg = uint8(lobe_bg) + uint8(applyTransf(double(hist_3D), ...
                double(lobe_mask(:,:,:,j)), params_2_refine(:,:,j), 'linear', 0));
        end   
    end
    if (size(lobe_ids,2)>1)
        lobe_bg = lobe_bg./(size(lobe_ids,2)-1);
    end

    [rigid_params(:,:,i), ~] = registerSlices (params_2_refine(:,:,i), hist_3D, ...
        lobe_mask(:,:,:,i), target, N, K, score_forward, lobe_bg);
end

fn = [path_out, 'params.mat'];
disp (['Saving params to ', fn])
save(fn , 'rigid_params');

hist_3D = hist_3D.*uint8(sum(mask_4D,4)>0);
hist_3D_R = hist_3D_R.*uint8(sum(mask_4D,4)>0);
hist_3D_G = hist_3D_G.*uint8(sum(mask_4D,4)>0);
hist_3D_B = hist_3D_B.*uint8(sum(mask_4D,4)>0);

%%%
%%% apply the transform and postprocess
%%%
reconstruction = applyTransf(hist_3D, lobe_mask>0, rigid_params, 'linear', 0);
reconstruction(reconstruction==0) = max(reconstruction(:));

reconstructionR = applyTransf(hist_3D_R, lobe_mask>0, rigid_params, 'linear', 0);
reconstructionR(reconstructionR==0) = max(reconstructionR(:));

reconstructionG = applyTransf(hist_3D_G, lobe_mask>0, rigid_params, 'linear', 0);
reconstructionG(reconstructionG==0) = max(reconstructionG(:));

reconstructionB = applyTransf(hist_3D_B, lobe_mask>0, rigid_params, 'linear', 0);
reconstructionB(reconstructionB==0) = max(reconstructionB(:));

mask_4D_reg = applyTransf(mask_4D, lobe_mask>0, rigid_params, 'nearest',0);


%%%
%%% Correct for Matlabs' weird orientation system
%%%
reconstruction  = postprocessRecontructionOutput(reconstruction);
reconstructionR = postprocessRecontructionOutput(reconstructionR);
reconstructionG = postprocessRecontructionOutput(reconstructionG);
reconstructionB = postprocessRecontructionOutput(reconstructionB);
mask_4D_reg     = postprocessRecontructionOutput(mask_4D_reg);


end
