%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstuct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear all;

this_script = mfilename('fullpath');
[this_path,name,ext] = fileparts(this_script);

%location of the util code
addpath(genpath([this_path,'\util\']))

%location of the data for the example provided here
path_base = [this_path, '\..\'];

for flipZ = 0:1
    % elastix parameter file
    param = [path_base,'\elastix\params\affine.txt'];

    % fixed image
    fIm = [path_base,'\example',example_id,'\output',output_id,'\step3_refineHistology\', ...
        num2str(flipZ),'\mha\lesion_masked.mha'];

    %moving image - crop volume to nodule and masked
    mIm = [path_base,'\example',example_id,'\data\CT\CT_crop_masked.mha'];

    %moving image - entire CT
    mIm2 = [path_base,'\example',example_id,'\data\CT\CT.mha'];

    %label of the nodule in the moving image - entire CT
    mIm3 = [path_base,'\example',example_id,'\data\CT\CT_label.mha'];

    Moving    = 'moving';
    Moving2   = 'uncropped';
    Mask      = 'mask';

    outParam  = 'TransformParameters.0.txt';

    %
    % Affine Registration
    %
    t0=[ path_base, '\example',example_id,'\output',output_id,'\step2_exhaustiveSearch\', num2str(flipZ),...
        '\1\', outParam ];

    path_out = [path_base, '\example',example_id,'\output',output_id,'\step4_refineCTToHistology\', ...
        num2str(flipZ),'\affine\']; 
    systemRunElastix(elastix_exe, fIm, mIm, t0, param, path_out);

    %apply same transform
    path_out_tr = [path_out, '\', Moving,'\'];
    systemRunTransformix(transformix_exe, mIm,[ path_out,'\', outParam],path_out_tr);
    mIm = [path_out_tr,'result.mha'];

    path_out_tr = [path_out, '\', Moving2,'\']; 
    systemRunTransformix(transformix_exe, mIm2,[ path_out,'\', outParam],path_out_tr);
    mIm2 = [path_out_tr,'result.mha'];

    path_out_tr = [path_out, '\', Mask,'\']; 
    systemRunTransformix(transformix_exe, mIm3,[ path_out,'\', outParam],path_out_tr);
    mIm3 = [path_out_tr,'result.mha'];



    %
    % Deformable
    %
    param = [path_base,'\elastix\params\deformable.txt'];

    t0='';
    path_out = [path_base, '\example',example_id,'\output',output_id,'\step4_refineCTToHistology\', ...
        num2str(flipZ),'\deformable\']; 
    systemRunElastix(elastix_exe, fIm, mIm, t0, param, path_out);

    %apply same transform
    path_out_tr = [path_out, '\', Moving,'\'];
    systemRunTransformix(transformix_exe, mIm,[ path_out,'\', outParam],path_out_tr);

    path_out_tr = [path_out, '\', Moving2,'\']; 
    systemRunTransformix(transformix_exe, mIm2,[ path_out,'\', outParam],path_out_tr);

    path_out_tr = [path_out, '\', Mask,'\']; 
    systemRunTransformix(transformix_exe, mIm3,[ path_out,'\', outParam],path_out_tr);

end
