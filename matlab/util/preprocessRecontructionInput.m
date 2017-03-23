function out_vol = preprocessRecontructionInput(in_vol, magni, rota, flip_LR_vol)

% resize images if they are still to large
if (magni~=1)
    out_vol = imresize(in_vol, 1/magni);
end

if (rota)
    if (numel(size(out_vol))==3)
        for (i=1:size(out_vol,3))
            out_vol(:,:,i) = rot90( rot90( out_vol(:,:,i)));
        end
    end
    
    if (numel(size(out_vol))==4)
        for (i=1:size(out_vol,4))
            out_vol(:,:,:,i) = rot90( rot90( out_vol(:,:,:,i)));
        end
    end
end

if (flip_LR_vol)
    disp 'Do flip';
    if (numel(size(out_vol))==3)
        for (i=1:size(out_vol,3))
            out_vol(:,:,i) = flipud( out_vol(:,:,i));
        end
    end
    
    if (numel(size(out_vol))==4)
        for (i=1:size(out_vol,4))
            out_vol(:,:,:,i) = flipud( out_vol(:,:,:,i));
        end
    end

end

end