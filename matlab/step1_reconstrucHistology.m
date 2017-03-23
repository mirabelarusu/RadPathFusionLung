%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstuct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear all;

this_script = mfilename('fullpath');
[this_path,name,ext] = fileparts(this_script);

%location of the util code
addpath(genpath([this_path,'\util\']))

%location of the data for the example provided here
path_base = [this_path, '\..\example',example_id,'\'];

% location downsample histology images
path_img = [path_base,'\data\Histology\imgs\']; 

% masks have to have the exact same name as the histology files
path_masks = [path_base,'\data\Histology\masks\']; 

% folder would be overwrote if it already exist
path_out = [path_base, '\output',output_id,'\step1_reconstructHistology\']; 

if (exist(path_out)==0)
    disp(['Creating folder', path_out]);
    mkdir(path_out)
else
    disp(['Writing to folder: "', path_out,'"']);
end

%how much to downsample
magni = 10;

[reconstruction, reconstructionR, reconstructionG, reconstructionB, mask_4D, rigid_params] = ...
    reconstruct ( path_img, path_out, 'path_masks', path_masks, ...
        'prefix_masks', prefix_masks, 'new_size', new_size, ...
        'magni', magni, 'image_order', image_order,'hist_prefix',hist_prefix,...
        'mask_prefix',mask_prefix, 'hist_extension', extension, 'id_roi', id_roi,...
        'id_sample', id_sample,'weight_entire_sample',weight_entire_sample);
    
    
     
% Final reconstruction voxel size: default 1,1,1: should reflect all
% scaling, magnification etc ...
page7 = 0.00372; %milimeters for page 7 (extracted)
shrink = 0.05; % percent shrinking of histology 
voxel_size = [page7*magni*(1+shrink),page7*magni*(1+shrink),spaceBetweenSlices];

paddingOnZ = 4;

% should the data be flipped on Z
flipOnZ = 0;

path_out_mha = [path_out, '/mha/', num2str(flipOnZ),'/'];
writeReconstruction(path_out_mha, prefix_masks,reconstruction, reconstructionR, ...
    reconstructionG, reconstructionB, mask_4D, voxel_size, flipOnZ, paddingOnZ, id_roi);


%should the data be flipped on Z
flipOnZ = 1;

path_out_mha = [path_out, '/mha/', num2str(flipOnZ),'/'];
writeReconstruction(path_out_mha, prefix_masks,reconstruction, reconstructionR, ...
    reconstructionG, reconstructionB, mask_4D, voxel_size, flipOnZ, paddingOnZ,id_roi);

