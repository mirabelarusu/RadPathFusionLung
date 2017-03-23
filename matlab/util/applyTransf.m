function out_vol = applyTransf(in_vol, in_masks, params, interpolator, bg)

out_vol = uint8(zeros(size(in_vol)));
denominator = uint8(zeros(size(in_vol)));

if (numel(size(in_vol))==3)
    for (i=1:size(in_vol,3))
        for (k = 1:size(params,3))
            trans_vol = uint8(imdef(params(i,:,k),double(in_vol(:,:,i)).*...
                double(in_masks(:,:,i,k)), interpolator,bg));
            
            out_vol(:,:,i) = out_vol(:,:,i) + trans_vol; 
            denominator(:,:,i) = denominator(:,:,i) + uint8(trans_vol>0);
        end
    end
end

if (numel(size(in_vol))==4)
    for (j =1:size(in_vol,3))
        for (i=1:size(in_vol,4))
            for (k = 1:size(params,3))
                trans_vol = uint8(imdef(params(j,:,k), double(in_vol(:,:,j,i)).*...
                    double(in_masks(:,:,j,k)), interpolator,bg));
                out_vol(:,:,j,i) = out_vol(:,:,j,i) + trans_vol;
                
                denominator(:,:,j,i) = denominator(:,:,j,i) + uint8(trans_vol>0);
            end
        end
    end
end

out_vol = out_vol./denominator;

end