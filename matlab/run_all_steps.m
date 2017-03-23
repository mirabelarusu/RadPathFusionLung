example_id = '';

%elastix executable; Dont forget \" to escape the space and allow elastix
%running
elastix_base = '"C:\Program Files\elastix_v4.8\\"';
elastix_exe = [elastix_base, 'elastix']; 
transformix_exe = [elastix_base, 'transformix']; 


% all important marks
prefix_masks = {'\sample\';'\lesion\';'\invasive\';'\blood\'};

%which mask should be used to drive the registration (used to crop the img)
id_roi = 2;
id_sample = 1;

%size that fit them all 
new_size = [6000,10500]; 

%how much to downsample
magni = 10;

%searches for histology slices called: Hist1.tiff, Hist2.tiff
hist_prefix = 'Hist';
mask_prefix = hist_prefix;

% data is read in alfabetic order, but that might not be what one desires
% default should be order = 1:4, but the images should also have the same
% numbers
image_order = 1:5;

extension = 'tif';

weight_entire_sample = 0.100;

for spaceBetweenSlices = 3:4
    output_id = ['_',num2str(spaceBetweenSlices)];

    try
        step1_reconstrucHistology
        step2_exhaustiveSearch
        step3_refineHistology
        step4_refineCTToHistology
        step5_mapHistologyOntoCT
    catch
        disp('Something went wrong!')
    end
end

