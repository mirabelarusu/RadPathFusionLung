%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Refine the reconstruction using CT as constraint
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear all;


this_script = mfilename('fullpath');
[this_path,name,ext] = fileparts(this_script);

%location of the util code
addpath(genpath([this_path,'\util\']))

%location of the data
path_base = [this_path, '\..\example',example_id,'\'];

% location downsample histology images
path_img = [path_base,'\data\Histology\imgs\']; 

% masks have to have the exact same name as the histology files
path_masks = [path_base,'\data\Histology\masks\']; 

% folder would be overwrote if it already exist
param_fn = [path_base, '\output',output_id,'\step1_reconstructHistology\params.mat'];
path_out = [path_base, '\output',output_id,'\step3_refineHistology\']; 


weight_entire_sample = 0.0001;  
     
% Final reconstruction voxel size: default 1,1,1: should reflect all
% scaling, magnification etc ...
page7 = 0.00372; %milimeters for page 7 (extracted)
shrink = 0.05; % percent shrinking of histology 
voxel_size = [page7*magni*(1+shrink),page7*magni*(1+shrink),spaceBetweenSlices];

%PaddingOnZ - how much padding should be added when experting the mha
paddingOnZ = 4;


% should the data be flipped on Z
flip_LR_vol = 0;
path_target = [path_base,'\output',output_id,'\step2_exhaustiveSearch\',num2str(flip_LR_vol),'\1\moving\result.mha'];
path_out_run = [path_out,'\',num2str(flip_LR_vol),'\']; 
if (exist(path_out_run)==0)
    disp(['Creating folder', path_out_run]);
    mkdir(path_out_run)
else
    disp(['Writing to folder: "', path_out_run,'"']);
end
load(param_fn)

[reconstruction, reconstructionR, reconstructionG, reconstructionB, mask_4D, rigidparams] = ...
    reconstruct ( path_img, path_out_run, 'path_masks', path_masks, ...
        'prefix_masks', prefix_masks, 'new_size', new_size, ...
        'hist_prefix',hist_prefix ,'mask_prefix', mask_prefix, 'hist_extension', extension, ....
        'id_roi', id_roi, 'id_sample', id_sample, 'weight_entire_sample',weight_entire_sample,...
        'path_target', path_target, 'ignored_padding_target', paddingOnZ,...
        'magni', magni, 'image_order', image_order, 'params', rigid_params);




path_out_mha = [path_out_run, '/mha/'];
writeReconstruction(path_out_mha, prefix_masks,reconstruction, reconstructionR, ...
    reconstructionG, reconstructionB, mask_4D, voxel_size, flip_LR_vol, paddingOnZ, id_roi);


load(param_fn)
flip_LR_vol = 1;
path_target = [path_base,'\output',output_id,'\step2_exhaustiveSearch\',num2str(flip_LR_vol),'\1\moving\result.mha'];
path_out_run = [path_out,'\',num2str(flip_LR_vol),'\']; 
if (exist(path_out_run)==0)
    disp(['Creating folder', path_out_run]);
    mkdir(path_out_run)
else
    disp(['Writing to folder: "', path_out_run,'"']);
end
[reconstruction, reconstructionR, reconstructionG, reconstructionB, mask_4D, rigidparams] = ...
    reconstruct ( path_img, path_out_run, 'path_masks', path_masks, ...
        'prefix_masks', prefix_masks, 'new_size', new_size, ...
        'hist_prefix',hist_prefix ,'mask_prefix', mask_prefix, 'hist_extension', extension, ....
        'id_roi', id_roi, 'flip_LR_vol', flip_LR_vol,...
        'id_sample', id_sample,'weight_entire_sample',weight_entire_sample,...
        'path_target', path_target, 'ignored_padding_target', paddingOnZ,...
        'magni', magni, 'image_order', image_order, 'params', rigid_params);


% should the data be flipped on Z
path_out_mha = [path_out_run, '/mha/'];
writeReconstruction(path_out_mha, prefix_masks,reconstruction, reconstructionR, ...
    reconstructionG, reconstructionB, mask_4D, voxel_size, flip_LR_vol, paddingOnZ, id_roi);

