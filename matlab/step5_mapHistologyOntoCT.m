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
outParam  = 'TransformParameters.0.txt';

for flipZ = 0:1
    %
    % Identify the inverse of the affine transform
    %
    param = [path_base,'\elastix\params\inverseAffine.txt'];

    % fixed image
    fIm = [path_base,'\example',example_id,'\data\CT\CT_crop_oversampled.mha'];

    %moving image - the reconstructed histology
    mIm = [path_base,'\example',example_id,'\output',output_id,'\step3_refineHistology\', num2str(flipZ),...
        '\mha\lesion_masked.mha'];

    t0=[ path_base, '\example',example_id,'\output',output_id,'\step4_refineCTToHistology\', ...
        num2str(flipZ),'\affine\', outParam ];

    path_out = [path_base, '\example',example_id,'\output',output_id,'\step5_mapHistologyOntoCT\', ...
        num2str(flipZ),'\inverse_affine\']; 
    systemRunElastix(elastix_exe, fIm, mIm, t0, param, path_out);

    %
    % Identify the inverse of the deformable transform
    %
    param = [path_base,'\elastix\params\inverseDeformable.txt'];

    % fixed image
    fIm = [path_base,'\example',example_id,'\output',output_id,'\step4_refineCTToHistology\', ...
        num2str(flipZ),'\affine\uncropped\result.mha'];

    t0=[ path_base, '\example',example_id,'\output',output_id,'\step4_refineCTToHistology\', ...
        num2str(flipZ),'\deformable\', outParam ];

    path_out = [path_base, '\example',example_id,'\output',output_id,'\step5_mapHistologyOntoCT\', ...
        num2str(flipZ),'\inverse_deformable\']; 
    systemRunElastix(elastix_exe, fIm, mIm, t0, param, path_out);

    %
    % Apply the same transform to the following 
    %
    mIm2 = [path_base,'\example',example_id,'\output',output_id,'\step3_refineHistology\', num2str(flipZ),...
        '\mha\reco_Gray.mha'];
    outPre2 = 'reconstruction_grayscale';
    path_out_tr2 = [path_out, '\', outPre2,'\'];
    systemRunTransformix(transformix_exe, mIm2,[ path_out,'\', outParam],path_out_tr2);

    mIm3 = [path_base,'\example',example_id,'\output',output_id,'\step3_refineHistology\', num2str(flipZ),...
        '\mha\blood_label.mha'];
    outPre3 = 'blood_label';
    path_out_tr3 = [path_out, '\', outPre3,'\'];
    systemRunTransformix(transformix_exe, mIm3,[ path_out,'\', outParam],path_out_tr3);

    mIm4 = [path_base,'\example',example_id,'\output',output_id,'\step3_refineHistology\', num2str(flipZ),...
        '\mha\invasive_label.mha'];
    outPre4 = 'invasion_label';
    path_out_tr4 = [path_out, '\', outPre4,'\'];
    systemRunTransformix(transformix_exe, mIm4,[ path_out,'\', outParam],path_out_tr4);

    mIm5 = [path_base,'\example',example_id,'\output',output_id,'\step3_refineHistology\', num2str(flipZ),...
        '\mha\lesion_masked.mha'];
    outPre5 = 'lesion';
    path_out_tr5 = [path_out, '\', outPre5,'\'];
    systemRunTransformix(transformix_exe, mIm5,[ path_out,'\', outParam],path_out_tr5);

    mIm6 = [path_base,'\example',example_id,'\output',output_id,'\step3_refineHistology\', num2str(flipZ),...
        '\mha\sample_label.mha'];
    outPre6 = 'sample_label';
    path_out_tr6 = [path_out, '\', outPre6,'\'];
    systemRunTransformix(transformix_exe, mIm6,[ path_out,'\', outParam],path_out_tr6);
    
    mIm7 = [path_base,'\example',example_id,'\output',output_id,'\step3_refineHistology\', num2str(flipZ),...
        '\mha\lesion_label.mha'];
    outPre7 = 'lesion_label';
    path_out_tr7 = [path_out, '\', outPre7,'\'];
    systemRunTransformix(transformix_exe, mIm7,[ path_out,'\', outParam],path_out_tr7);


    % once the inverse of the affine is identified, it needs to be edited
    % to allow proper applying it. Here we removed the parameter: 
    % InitialTransformParametersFileName with "NoInitialTransform"
    path_in = [path_base, '\example',example_id,'\output',output_id,'\step5_mapHistologyOntoCT\', ...
        num2str(flipZ),'\inverse_affine\']; 

    
    fn_param_in = [path_in,outParam]
    fid_in = fopen(fn_param_in)

    % get path and extension
    [pathstr,name,ext] = fileparts(fn_param_in);    
    fn_param_out = [pathstr,'\',name,'_modified',ext];
    fid_out = fopen(fn_param_out,'w+t');

    tline = fgetl(fid_in);
    lineNofile = 0;
    while ischar(tline)
        if (lineNofile == 3)
            tline = '(InitialTransformParametersFileName "NoInitialTransform")';
        end
        fprintf(fid_out,'%s\n',tline);

        lineNofile = lineNofile+1;
        tline = fgetl(fid_in);
    end

    fclose(fid_in);
    fclose(fid_out);


    %
    % apply the inverse affine transform to the already inversed derformation
    %
    path_out = [path_base, '\example',example_id,'\output',output_id,'\step5_mapHistologyOntoCT\', ...
        num2str(flipZ),'\inverse_deformable_and_affine\']; 

    mIm2 = [path_out_tr2,'\result.mha'];
    path_out_tr = [path_out, '\', outPre2,'\'];
    systemRunTransformix(transformix_exe, mIm2,fn_param_out,path_out_tr);
    try
        copyfile([path_out_tr, '\result.mha'], [path_out_tr, '..\..\',outPre2,'.mha'])
    catch 
        disp(['Something went wrong when copying file', path_out_tr, '\result.mha'] )
    end

    mIm3 = [path_out_tr3,'\result.mha'];
    path_out_tr = [path_out, '\', outPre3,'\'];
    systemRunTransformix(transformix_exe, mIm3,fn_param_out,path_out_tr);
    try
        copyfile([path_out_tr, '\result.mha'], [path_out_tr, '..\..\',outPre3,'.mha'])
    catch 
        disp(['Something went wrong when copying file', path_out_tr, '\result.mha'] )
    end

    mIm4 = [path_out_tr4,'\result.mha'];
    path_out_tr = [path_out, '\', outPre4,'\'];
    systemRunTransformix(transformix_exe, mIm4,fn_param_out,path_out_tr);
    try
        copyfile([path_out_tr, '\result.mha'], [path_out_tr, '..\..\',outPre4,'.mha'])
    catch 
        disp(['Something went wrong when copying file', path_out_tr, '\result.mha'] )
    end

    mIm5 = [path_out_tr5,'\result.mha'];
    path_out_tr = [path_out, '\', outPre5,'\'];
    systemRunTransformix(transformix_exe, mIm5,fn_param_out,path_out_tr);
    try
        copyfile([path_out_tr, '\result.mha'], [path_out_tr, '..\..\',outPre5,'.mha'])
    catch 
        disp(['Something went wrong when copying file', path_out_tr, '\result.mha'] )
    end

    mIm6 = [path_out_tr6,'\result.mha'];
    path_out_tr = [path_out, '\', outPre6,'\'];
    systemRunTransformix(transformix_exe, mIm6,fn_param_out,path_out_tr);
    try
        copyfile([path_out_tr, '\result.mha'], [path_out_tr, '..\..\',outPre6,'.mha'])
    catch 
        disp(['Something went wrong when copying file', path_out_tr, '\result.mha'] )
    end
    
    mIm7 = [path_out_tr7,'\result.mha'];
    path_out_tr = [path_out, '\', outPre7,'\'];
    systemRunTransformix(transformix_exe, mIm7, fn_param_out,path_out_tr);
    try
        copyfile([path_out_tr, '\result.mha'], [path_out_tr, '..\..\',outPre7,'.mha'])
    catch 
        disp(['Something went wrong when copying file', path_out_tr, '\result.mha'] )
    end
end