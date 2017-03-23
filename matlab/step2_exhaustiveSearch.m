%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Full search to identify the position of CT relative
% to histology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear all;

this_script = mfilename('fullpath');
[this_path,name,ext] = fileparts(this_script);

%location of the util code
addpath(genpath([this_path,'\util\']))



%location of the data

%location of the data
path_base = [this_path, '\..\'];


for flipZ = 0:1
    % elastix parameter file
    param = [path_base,'\elastix\params\fullSearch.txt'];

    % fixed image
    fIm = [path_base,'\example',example_id,'\output',output_id,'\step1_reconstructHistology\mha\', ...
        num2str(flipZ),'\lesion_masked.mha'];

    %moving image - crop volume to nodule and masked
    mIm = [path_base,'\example',example_id,'\data\CT\CT_crop_masked.mha'];

    %moving image - entire CT
    mIm2 = [path_base,'\example',example_id,'\data\CT\CT.mha'];

    Mask      = 'moving';
    Moving    = 'uncropped';
    outParam  = 'TransformParameters.0.txt';

    %
    % Exhaustive Search
    %
    path_out = [path_base, '\example',example_id,'\output',output_id,'\step2_exhaustiveSearch\', ...
        num2str(flipZ),'\0\']; 
    systemRunElastix(elastix_exe, fIm, mIm, '', param, path_out);

    %apply same transform
    path_out_tr = [path_out, '\', Mask,'\'];
    systemRunTransformix(transformix_exe, mIm,[ path_out,'\', outParam],path_out_tr);

    path_out_tr = [path_out, '\', Moving,'\']; 
    systemRunTransformix(transformix_exe, mIm2,[ path_out,'\', outParam],path_out_tr);


    %
    % Affine Refinement of the found solutions
    %
    % elastix parameter file; different than before
    param = [path_base,'\elastix\params\affinePostFullSearch.txt'];
    t0=[path_out,'\', outParam ];
    path_out = [path_base, '\example',example_id,'\output',output_id,'\step2_exhaustiveSearch\', ...
        num2str(flipZ),'\1\']; 
    systemRunElastix(elastix_exe, fIm, mIm, t0, param, path_out);

    %apply same transform
    path_out_tr = [path_out, '\', Mask,'\'];
    systemRunTransformix(transformix_exe, mIm,[ path_out,'\', outParam],path_out_tr);

    path_out_tr = [path_out, '\', Moving,'\']; 
    systemRunTransformix(transformix_exe, mIm2,[ path_out,'\', outParam],path_out_tr);
end

