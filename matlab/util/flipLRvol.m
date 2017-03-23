function flippedVol = flipLRVol(vol) 

s = size(vol);

if (numel( s )==3)
    flippedVol = zeros(s);
    for i = 1:s(3)
        flippedVol(:,:,i) = fliplr(vol(:,:,i));
    end
else
    disp('Input must be a 3D volume!')
    flippedVol = vol;
end

end