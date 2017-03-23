function [rigid_params, d_reg] = registerSlices(params_2_refine, vol, mask, ...
    target, Nbins, Kslices, score_forward, bg)

if (size(vol,3)==1)
    disp(['There must be at least ', num2str(Kslices), 'in the volume to register'])
    return;
end

res1 = fspecial('gaussian',[50 50],20);
res2 = fspecial('gaussian',[20 20],10);
d_reg = uint8(zeros(size(vol)));
masked_vol = uint8(double(vol).*double(mask));

if (size(params_2_refine,1)==0)
    rigid_params = zeros(size(vol,3),5);
    rigid_params(1,4) = 1;
    rigid_params(1,5) = 1;
else
    rigid_params = params_2_refine;
end;

d_reg  = masked_vol;

for (i = 1:size(vol,3))

    if (sumall(masked_vol(:,:,i))==0) % is empty
        disp(['Noting to register at ', num2str(i)])
        continue
    end
    
    if (i==1)
        weightRandom=0;
    else
        weightRandom=1;
    end
        
    if (mod(i-1,1)==0)
        noDisp=false;
    else
        noDisp=true;
    end
    

   t = createReconstructionTarget(d_reg, i, target, Kslices, score_forward,rigid_params);
    
    clear t1;
    clear d_reg1;  
    for (idxT=1:size(t,3))
        t1(:,:,idxT) = imresize(imfilter(t(:,:,idxT), res1, 'same'),1.0/4.0);
    end
     
    if (size(t,3)>Kslices)
        %removed the mri from the estimation
        stripped_t = t;
        stripped_t(:,:,Kslices+1) = 0;
    else
        stripped_t = t;
    end
        
    % allows to compute the Center of mass of t
        
    t_proj_mask = sum(stripped_t,3)>0;
    t_mask_ch = bwconvhull(t_proj_mask);
    t_rp = regionprops(t_mask_ch,{'Orientation','Centroid'});
    com_t = cat(1,t_rp.Centroid);
    orient_t = cat(1,t_rp.Orientation);
    
    if (size(orient_t,1)>0) % a target exist    
        rigid_params(i,4) = 1;
        rigid_params(i,5) = 1;
        
        % get the center of mass of the transformed image
        
        moved_mask = masked_vol(:,:,i)>0;
        m_mask_ch = bwconvhull(moved_mask);
        m_rp = regionprops(m_mask_ch,{'Orientation','Centroid'});
        com_m = cat(1,m_rp.Centroid);
        orient_m = cat(1,m_rp.Orientation);
        
        % apply rotation only which will give a better estimate of the
        % translation
        ang_diff = orient_t-orient_m;
        if (ang_diff<-90) % aligns the dominant axes the wrong way
            ang_diff = 180 + ang_diff;
        end
        
        if (ang_diff>90) % aligns the dominant axes the wrong way
            ang_diff = 180 - ang_diff;
        end
        newparam(1) = (ang_diff*pi/180.0);
        newparam(2) = 0;
        newparam(3) = 0;
        newparam(4) = 1;
        newparam(5) = 1;
        moved_mask = imdef(newparam, masked_vol(:,:,i), 'linear', 0)>0;
        m_mask_ch = bwconvhull(moved_mask);
        m_rp = regionprops(m_mask_ch,{'Orientation','Centroid'});
        com_m = cat(1,m_rp.Centroid);
        
        
        %%%
        %%% First layer in the pyramid
        %%%
        newparam(1) = (ang_diff/180.0*pi+rigid_params(i,1))/2.0;
        newparam(2) = ((com_t(1)-com_m(1))+rigid_params(i,2))/8.0;
        newparam(3) = ((com_t(2)-com_m(2))+rigid_params(i,3))/8.0;
        %newparam(4) = 1;
        %newparam(5) = 1;
                
        %imshow([moved_mask+moved_mask2+t_proj_mask moved_mask2+t_proj_mask])
        for (idx =1:3)
            newoffset(idx) = newparam(idx);
        end
        
        
        moved_mask = imdef(newparam, masked_vol(:,:,i), 'linear', 0)>0; 
        
        %if (maxall(moved_mask+t_proj_mask)==1) % they don't overlap; so just use the old params 
            disp(['Something went wrong when estimating rough alignment'])
            newparam(1) = rigid_params(i,1);
            newparam(2) = rigid_params(i,2)/4.0;
            newparam(3) = rigid_params(i,3)/4.0;

            for (idx =1:3)
                newoffset(idx) = newparam(idx);
            end
        %end            

    else
        newparam(1) = rigid_params(i,1);
        newparam(2) = rigid_params(i,2)/4.0;
        newparam(3) = rigid_params(i,3)/4.0;
        
        for (idx =1:3)
            newoffset(idx) = newparam(idx);
        end
    end
    
    disp(rigid_params)

    m = 1/4.0; % downsample 4 times
    fixed = uint32(t1);
    moving = uint32(imresize(imfilter(masked_vol(:,:,i), res1, 'same'),m));
    moving_bg = uint32(imresize(imfilter(bg(:,:,i), res1, 'same'),m));
    [d_reg1(:,:,i), rigidparam] =  rigidreg2_groupwise(fixed, moving, ...
        Nbins, 'NoDisp', noDisp, 'Measure','MMI', 'scales', 0,...  
        'paramoffsets', newoffset, 'Itarget_bg', moving_bg, 'number_many_slices',4);
    
    clear t1;
    clear d_reg1;
    for (idxT=1:size(t,3))
        t1(:,:,idxT) = imresize(imfilter(t(:,:,idxT), res2, 'same'),1.0/2.0);
    end
    
    %%%
    %%% Second layer in the pyramid
    %%%
    newparam(1) = rigidparam(1);
    newparam(2) = rigidparam(2)*2;
    newparam(3) = rigidparam(3)*2;
    
    for (idx =1:3)
        newoffset(idx) = newparam(idx);
    end

    
    m = 1/2.0; % downsample 2 times
    fixed = uint32(t1);
    moving = uint32(imresize(imfilter(masked_vol(:,:,i), res2, 'same'),m));
    moving_bg = uint32(imresize(imfilter(bg(:,:,i), res2, 'same'),m));    
    [d_reg1(:,:,i), rigidparam] =  rigidreg2_groupwise(fixed, moving,...
        Nbins, 'NoDisp', noDisp, 'Measure','MMI', 'scales', 0,... 
        'paramoffsets', newoffset, 'Itarget_bg', moving_bg, 'number_many_slices',4);

    %%%
    %%% Third layer in the pyramid
    %%%
    newparam(1) = rigidparam(1);
    newparam(2) = rigidparam(2)*2;
    newparam(3) = rigidparam(3)*2;
    
    for (idx =1:3)
        newoffset(idx) = newparam(idx);
    end
    
    fixed = uint32(t);
    moving = uint32(masked_vol(:,:,i));
    moving_bg = bg(:,:,i);    
    [d_reg(:,:,i), rigidparam] = rigidreg2_groupwise(fixed, moving, ...
        Nbins, 'NoDisp', noDisp, 'Measure','MMI', 'scales', 0, ...
        'paramoffsets', newoffset, 'Itarget_bg', moving_bg, 'number_many_slices',4);
    
    rigid_params(i,:) = rigidparam;

end

end
