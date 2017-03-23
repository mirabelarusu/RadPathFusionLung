function [hist_3D, hist_3D_R, hist_3D_G, hist_3D_B, mask_4D, target] = ...
    readReconstructInput (path_hist, hist_prefix, mask_prefix, image_order, ...
    new_size, path_masks, prefix_masks, path_target, N, do_threshold, ...
    padding, extension, flip_LR_vol)
% readSlices - reads the input, histology including annotations if existant
% do_threshold - should use otsu threshold on histology to remove the grey
% background?
% MR
files = dir(path_hist);
%numberImages = size(files,1);
numberImages = numel(image_order);


%if order was not provided, use the order as 
if (size(image_order)==0)
    image_order = 1:(size(files,1));
end
%image_order = [1,2, image_order+2];

hist_3D = uint8(zeros(new_size(1),new_size(2),numberImages));
hist_3D_R = hist_3D;
hist_3D_G = hist_3D;
hist_3D_B = hist_3D;

mask_4D = uint8(zeros(new_size(1),new_size(2),numberImages,size(prefix_masks,1)));

for (x = 1:numberImages)
    if (x > numel(image_order))
        disp(['Skipping the rest of images >', num2str(x)]);
        break
    end
    i = image_order(x);
    
    %histology intensity
    fn = [path_hist, hist_prefix, num2str(i),'.', extension];
    
    %read histology
    hist_2D_color   = imread(fn);
    disp(['Done opening ', int2str(i) ,' ', fn, ' at [',...
        num2str(size(hist_2D_color,1)),',',num2str(size(hist_2D_color,2)),']']);
    
    if (size(hist_2D_color,3)>2)
        hist_2D         = rgb2gray(hist_2D_color(:,:,1:3));
    else
        hist_2D         = rgb2gray(hist_2D_color);
    end;
    
    if (new_size(1)>0)
        hist_2D         = impad(hist_2D, new_size,cast(1,class(hist_2D)));
        hist_2D_color   = impad(hist_2D_color, new_size,cast(1,class(hist_2D)));
    end
    
%     if (do_threshold)% consider the entire sample and then threshold
%         thresh = graythresh(hist_2D);
%         disp(['Found threshold ', num2str(thresh)])
%         for (j = 1:3)
%             hh = hist_2D_color(:,:,j);
%             hh(hist_2D>thresh) = 1;
%             hist_2D_color(:,:,j) = hh;
%         end
%         hist_2D(hist_2D>thresh)=1;
%     else
%         thresh = 0;
%     end
    
    hist_2D = rescale(hist_2D)*(N-1);
    hist_2D_color = rescale(hist_2D_color)*(N-1);
    
    hist_3D(:,:,x) = hist_2D;
    hist_3D_R(:,:,x) = hist_2D_color(:,:,1);
    hist_3D_G(:,:,x) = hist_2D_color(:,:,2);
    hist_3D_B(:,:,x) = hist_2D_color(:,:,3);
    
    
    
    
    %read all the masks
    for (j = 1:size(prefix_masks,1))
        fn = [path_masks, char(prefix_masks(j)), '/', mask_prefix, num2str(i),'.', extension];
        disp(['  Mask_', int2str(j) ,' ', fn]);
        
        try
            mask = imread(fn);
            
            if (size(mask,3)>2)
                mask = rgb2gray(mask(:,:,1:3));
            end;
            
            if (new_size(1)>0)
                mask        = impad(mask, new_size, cast(0,class(mask)));
            end
        catch
            disp(['Can not load', fn,'!!! Setting it to zero!'])
            mask = zeros(new_size);
        end
        mask_4D(:,:,x,j) = mask;
    end

end

% read target if existant: 
target = 0;  
if (size(path_target,1) > 0 && exist(path_target)) 
	target = mha_read_volume(path_target);
	target = rescale(target).*(N-1);
    new_target = target(:,:,(padding+1):(size(target,3)-padding));
    if (flip_LR_vol)   
        disp 'Do flip target';
        for (i=1:size(new_target,3))
            new_target(:,:,i) = fliplr(new_target(:,:,i));
        end
    end
    %if (image_order(1)>image_order(size(image_order,2)))
    %    target = new_target(:,:,size(new_target,3):-1:1);
    %else
        target = new_target;
    %end
    clear new_target
    disp(['Reading target path from "',path_target,'"  @', num2str(size(target))])
else
    disp(['Could not read "',path_target,'"! Check if the file exists. No target volume will be considered.'])
end

end