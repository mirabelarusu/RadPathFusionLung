function out_vol = postprocessRecontructionOutput(in_vol)

out_vol = zeros(size(in_vol,2),size(in_vol,1),size(in_vol,3),size(in_vol,4));

if (numel(size(in_vol))==3)
    for (i=1:size(in_vol,3))
        out_vol(:,:,i) = in_vol(:,:,i)';
    end
end

if (numel(size(in_vol))==4)
    for (j=1:size(in_vol,4))
        for (i=1:size(in_vol,3))
            out_vol(:,:,i,j) = in_vol(:,:,i,j)';
        end
    end
end
end