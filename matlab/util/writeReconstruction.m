function status = writeReconstruction(path_out, prefix_masks, reconstruction, reconstructionR, reconstructionG,...
    reconstructionB, mask_3D, voxel_size, flipOnZ, paddingOnZ, idCreateMask)

    if (nargin<=10)
        idCreateMask = 1;
    end

    % create folder for output if it doesn't exists
    if (exist(path_out)==0)
        disp(['Creating folder', path_out]);
        mkdir(path_out)
    else
        disp(['Writing to folder: "', path_out,'"']);
    end

    newVol = zeros(size(reconstruction,1),size(reconstruction,2), size(reconstruction,3)+2*paddingOnZ);
    
    newRange = (1+paddingOnZ):(size(reconstruction,3)+paddingOnZ);
    
    newVol(:,:,newRange) = reconstruction;
    if (flipOnZ==1) newVol = flipLRvol(newVol); end
    mha_write_volume([path_out, 'reco_Gray.mha'],  newVol, voxel_size);
    newVol_hist = newVol;
    
    for (j = 1:size(prefix_masks,1))
        newVol = zeros(size(mask_3D(:,:,:,j),1),size(mask_3D(:,:,:,j),2), size(mask_3D(:,:,:,j),3)+2*paddingOnZ);
        newVol(:,:,newRange) = mask_3D(:,:,:,j);
        if (flipOnZ==1)
           newVol(:,:,newRange) =  flipLRvol(mask_3D(:,:,:,j));
        end
        stri = char(prefix_masks(j));
        stri = stri(1:(length(stri)-1));
        mha_write_volume([path_out, stri, '_label.mha'], newVol, voxel_size);
        if (j == idCreateMask)
            mha_write_volume([path_out, stri, '_masked.mha'], (maxall(newVol_hist)-newVol_hist).*newVol, voxel_size);
        end
    end

status = 1;